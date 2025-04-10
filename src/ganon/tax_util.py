import time
import pandas as pd
import re
import os
import gzip

from ganon.util import download, run, print_log, rm_files
from io import StringIO


def parse_sequence_accession(input_files, info_cols, build_output_folder):
    """
    Parse sequence accession (anything from > to the first space) in all input files
    Split files by sequence and write them in subfolders
    """
    info = pd.DataFrame(columns=info_cols)

    # Create subfolders to distribute files
    n_folders = 10
    run("seq 0 {n_folders} | xargs -i mkdir -p \"{build_output_folder}/{{}}\"".format(n_folders=n_folders-1, build_output_folder=build_output_folder), shell=True)
    
    # Randomly distribute sequences into subfolders
    # cat | zcat -> compability with osx
    for file in input_files:
        run_cat = """
        cat {file} {zcat} | 
        awk '/^>/{{
                id=substr($1,2);
                subf=int({n_folders} * rand())"/";
                file=(\"{build_output_folder}\" subf id \".fna\"); 
                print id\"\t\"file;
            }}
            {{print $0 > file}}'
        """.format(file=file, 
                   zcat="| zcat" if file.endswith(".gz") else "", 
                   build_output_folder=build_output_folder, 
                   n_folders=n_folders)
        
        stdout = run(run_cat, shell=True, ret_stdout=True)
        tmp_info = pd.read_csv(StringIO(stdout), header=None, names=['target', 'file'], delimiter="\t")
        info = pd.concat([info, tmp_info])

    return info

def parse_file_accession(input_files, info_cols):
    """
    Look for genbank/refseq assembly accession* pattern in the filename
    if not found, return basename of the file as target

    *https://support.nlm.nih.gov/knowledgebase/article/KA-03451/en-us
    *https://https.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/gca-and-gcf-explained/
    """
    info = pd.DataFrame(columns=info_cols)

    assembly_accessions = []
    assembly_accession_pattern = re.compile("GC[A|F]_[0-9]+\.[0-9]+")
    for file in input_files:
        match = assembly_accession_pattern.search(file)
        assembly_accessions.append((match.group() if match else os.path.basename(file), file))
    info[["target", "file"]] = pd.DataFrame(assembly_accessions)

    return info

def parse_genome_size_files(cfg, build_output_folder):
    """
    Download and parse auxiliary files to determine approximate genome size for each taxa
    NCBI based on species_genome_size.txt and GTDB on _metadata.tar.gz files
    Returns sizes for leaf nodes in the provided taxonomy in a dict {node:size}
    If node has no size information, returns {node:0}
    """

    # Download or set files
    tx = time.time()
    if not cfg.genome_size_files:
        print_log("Downloading and parsing auxiliary files for genome size estimation", cfg.quiet)
        if cfg.taxonomy == "ncbi":
            files = download([cfg.ncbi_url + "/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz"], build_output_folder)
        elif cfg.taxonomy == "gtdb":
            files = download([cfg.gtdb_url + "/ar53_metadata.tsv.gz", cfg.gtdb_url + "/bac120_metadata.tsv.gz"], build_output_folder)
    else:
        print_log("Parsing auxiliary files for genome size", cfg.quiet)
        files = cfg.genome_size_files

    leaves_sizes = {}
    if cfg.taxonomy == "ncbi":
        for file in files:
            with gzip.open(file, "rt") as f:
                # skip first line wiht header
                # #species_taxid  min_ungapped_length  max_ungapped_length  expected_ungapped_length  number_of_genomes  method_determined
                next(f)
                for line in f:
                    fields = line.rstrip().split("\t")
                    leaves_sizes[fields[0]] = int(fields[3])

    elif cfg.taxonomy == "gtdb":
        for file in files:
            with gzip.open(file, "rt") as f:
                # skip first line wiht header
                # col 0: accession (with GC_ RF_ prefix), col 16: genome_size, col 19: gtdb_taxonomy (d__Archaea;p__Thermoproteota;...)
                next(f)
                for line in f:
                    fields = line.rstrip().split("\t")
                    t = fields[19].split(";")[-1]  # species taxid (leaf)
                    # In GTDB, several genome sizes are available for each node
                    # accumulate them in a list and make average
                    if t not in leaves_sizes:
                        leaves_sizes[t] = []
                    leaves_sizes[t].append(int(fields[16]))
                    
        # Average sizes
        for t in list(leaves_sizes.keys()):
            leaves_sizes[t] = int(sum(leaves_sizes[t])/len(leaves_sizes[t]))
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return leaves_sizes


def parse_genome_size_tax(tax_files):
    """
    Parse last column of a .tax file and retrieve genome_sizes to a dict {node: size}
    If nodes are repeated among files with different sizes, keep largest
    """
    genome_sizes = {}
    for f in tax_files:
        with open(f, "r") as file:
            for line in file:
                node, _, _, _, gsize = line.rstrip().split("\t")
                # keep largest genome size
                gsize = int(gsize)
                if node in genome_sizes and genome_sizes[node] > gsize:
                    continue
                genome_sizes[node] = gsize
    return genome_sizes


def get_genome_size(cfg, nodes, tax, build_output_folder):
    """
    Estimate genome sizes based on auxiliary files
    Only used nodes and lineage are calculated, based on the full set of values provided
    If information of a certain node is not provided, uses the closest estimate of parent nodes
    """
    genome_sizes = {}
    if cfg.skip_genome_size:
        # Skipping genome sizes, all set to 1
        for node in nodes:
            for t in tax.lineage(node):
                genome_sizes[t] = 1
    else:
        # Download and parse auxiliary files containing genome sizes
        leaves_sizes = parse_genome_size_files(cfg, build_output_folder)

        tx = time.time()
        print_log("Estimating genome sizes", cfg.quiet)

        # Check if entries are on tax and distribute values to available tax. leaves
        for t in list(leaves_sizes.keys()):
            if not tax.latest(t):
                del leaves_sizes[t]
            else:
                # Store genome size estimation for all leaf nodes available in the taxonomy
                for leaf in tax.leaves(t):
                    leaves_sizes[leaf] = leaves_sizes[t]

        # Calculate genome size estimates for used nodes (and their lineage)
        # using the complete content of leaves_sizes (keeping approx. the same estimates between different dbs)
        
        for node in nodes:
            # For the lineage of each target node
            for t in tax.lineage(node):
                # Skip if already calculated
                if t not in genome_sizes:
                    cnt = 0
                    avg = 0
                    # Make average of available genome sizes in children leaves
                    for leaf in tax.leaves(t):
                        if leaf in leaves_sizes:
                            cnt += 1
                            avg += leaves_sizes[leaf]
                    genome_sizes[t] = int(avg / cnt) if cnt else 0

        # If there is no matching between taxonomy and leaves, average the whole and save to root to be redistributed in the next step
        if sum(genome_sizes.values())==0:
            if leaves_sizes:
                genome_sizes[tax.root_node] = int(sum(leaves_sizes.values())/len(leaves_sizes))
            else:
                genome_sizes[tax.root_node] = 1
        # Check nodes without genome size info (0) and use closest value from parent lineage
        for node in nodes:
            if genome_sizes[node] == 0:
                # Fill lineage of zeros with latest genome size estimation
                for t in tax.lineage(node):
                    if genome_sizes[t] == 0:
                        genome_sizes[t] = genome_sizes[tax.parent(t)]

        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return genome_sizes


def get_file_info(cfg, info, tax, build_output_folder):
    if cfg.taxonomy == "ncbi" or (cfg.taxonomy == "skip" and cfg.level == "assembly"):
        assembly_summary_urls = []
        assembly_summary_files = []

        for assembly_summary in cfg.ncbi_file_info:
            # Either a file or prefix for download
            if assembly_summary in cfg.choices_ncbi_file_info:
                # split if _historical
                assembly_summary_urls.append(cfg.ncbi_url +
                                             "/genomes/" + assembly_summary.split("_")[0] + 
                                             "/assembly_summary_" + assembly_summary + ".txt")
            else:
                assembly_summary_files.append(assembly_summary)

        if assembly_summary_urls:
            tx = time.time()
            print_log("Downloading assembly_summary files", cfg.quiet)
            assembly_summary_files.extend(download(assembly_summary_urls, build_output_folder))
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        tx = time.time()
        print_log("Parsing assembly_summary files", cfg.quiet)
        count_assembly_summary = parse_assembly_summary(info, assembly_summary_files, cfg.level)
        for assembly_summary_file, cnt in count_assembly_summary.items():
            print_log(" - " + str(cnt) + " entries found in the " + assembly_summary_file.split("/")[-1] + " file", cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    elif cfg.taxonomy == "gtdb":
        tx = time.time()
        print_log("Parsing gtdb files", cfg.quiet)
        # update nodes to info
        info.update(get_gtdb_target_node(tax, cfg.level))
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)


def get_gtdb_target_node(tax, level):
    """
    Parse GTDB taxonomy files and return target (accessions) node (gtdb taxid)
    use taxonomy files from multitax tax.sources (set with output_prefix on GtdbTx)
    """
    gtdb_target_node = pd.DataFrame(columns=["node"])
    for source in tax.sources:
        # Parse trimming prefix RS_ or GB_
        gtdb_target_node = pd.concat([gtdb_target_node,
                                      pd.read_csv(source,
                                                  sep="\t",
                                                  names=["target", "node"],
                                                  header=None,
                                                  index_col="target",
                                                  converters={"target": lambda x: x[3:],
                                                              "node": lambda l: l.split(";")[-1]}
                                                  )
                                      ])

    # Create specialization if requested for --level
    if level == "assembly":
        gtdb_target_node["specialization"] = gtdb_target_node.index
        gtdb_target_node["specialization_name"] = gtdb_target_node["node"].apply(tax.name)

    return gtdb_target_node


def get_sequence_info(cfg, info, tax, build_output_folder):
    if cfg.taxonomy == "ncbi":
        # Max. sequences to use eutils in auto mode
        max_seqs_eutils = 50000

        if not cfg.ncbi_sequence_info:
            # auto
            if info.shape[0] > max_seqs_eutils:
                mode = ["nucl_gb", "nucl_wgs"]
            else:
                mode = ["eutils"]
        elif "eutils" in cfg.ncbi_sequence_info:
            mode = ["eutils"]
        else:
            mode = cfg.ncbi_sequence_info  # custom accession2taxid prefixes or files

        if mode[0] == "eutils":
            tx = time.time()
            print_log("Retrieving sequence information from NCBI e-utils", cfg.quiet)
            info.update(run_eutils(cfg, info, build_output_folder, skip_taxid=False, level=cfg.level))
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
        else:
            if tax:
                acc2txid_urls = []
                acc2txid_files = []
                for acc2txid in mode:
                    # Either a file or prefix for download
                    if acc2txid in cfg.choices_ncbi_sequence_info:
                        acc2txid_urls.append(cfg.ncbi_url + "/pub/taxonomy/accession2taxid/" +
                                             acc2txid + ".accession2taxid.gz")
                    else:
                        acc2txid_files.append(acc2txid)

                if acc2txid_urls:
                    tx = time.time()
                    print_log("Downloading accession2taxid files", cfg.quiet)
                    acc2txid_files.extend(download(acc2txid_urls, build_output_folder))
                    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

                tx = time.time()
                print_log("Parsing accession2taxid files", cfg.quiet)
                count_acc2txid = parse_acc2txid(info, acc2txid_files)
                for acc2txid_file, cnt in count_acc2txid.items():
                    print_log(" - " + str(cnt) + " entries found in the " +
                              acc2txid_file.split("/")[-1] + " file", cfg.quiet)
                print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

            # If assembly is requested, need to use e-utils
            if cfg.level == "assembly":
                print_log("Retrieving assembly/name information from NCBI e-utils", cfg.quiet)
                info.update(run_eutils(cfg, info, build_output_folder, skip_taxid=True, level="assembly"))
                print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    elif cfg.taxonomy == "gtdb":
        tx = time.time()
        # The only way to link gtdb to sequences is through NCBI assembly accessions
        print_log("Retrieving assembly/name information from NCBI e-utils", cfg.quiet)
        # Get NCBI assembly accession for every sequence
        assembly_info = run_eutils(cfg, info, build_output_folder, skip_taxid=True, level="assembly")
        # Get map of accessions from from gtdb taxonomy
        gtdb_target_node = get_gtdb_target_node(tax, cfg.level)
        # Update info with merged info
        info.update(assembly_info.join(on="specialization", other=gtdb_target_node, lsuffix="_acc"))
        # if no specialization was requested, delete data
        if cfg.level not in cfg.choices_level:
            info["specialization"] = None
            info["specialization_name"] = None
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    elif cfg.taxonomy == "skip":
        if cfg.level == "assembly":
            tx = time.time()
            print_log("Retrieving sequence information from NCBI e-utils", cfg.quiet)
            info.update(run_eutils(cfg, info, build_output_folder, skip_taxid=True, level="assembly"))
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)


def parse_acc2txid(info, acc2txid_files):
    count_acc2txid = {}
    unique_acc = set(info.index)

    # Parse very large acc2txid files in chunks to limit memory usage
    chunksize = 10 ** 6
    for acc2txid in acc2txid_files:
        count_acc2txid[acc2txid] = 0
        with pd.read_csv(acc2txid,
                         sep="\t",
                         header=None,
                         skiprows=1,
                         usecols=[1, 2],
                         names=["target", "node"],
                         index_col="target",
                         converters={"target": lambda x: x if x in unique_acc else None, "node": str},
                         chunksize=chunksize) as reader:

            for tmp_acc_node in reader:
                tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()]  # keep only seqids used
                tmp_acc_node = tmp_acc_node[tmp_acc_node["node"] != "0"]  # filter out taxid==0

                # If there was any match
                if tmp_acc_node.shape[0]:
                    # merge node(taxid) retrieved based on target(accesion)
                    info.update(tmp_acc_node)

                    # save count to return
                    count_acc2txid[acc2txid] += tmp_acc_node.shape[0]
                    
                    #if already found all seqids no need to parse all files till the end)
                    if sum(count_acc2txid.values()) == len(unique_acc):
                        break

    return count_acc2txid


def parse_assembly_summary(info, assembly_summary_files, level):
    count_assembly_summary = {}
    unique_acc = set(info.index)

    for assembly_summary in assembly_summary_files:
        # Detect header lines beforehand, so pandas.read_csv can read it properly
        header_lines = 0
        with open(assembly_summary, 'r') as ass_sum:
            for line in ass_sum:
                if line[0] == "#":
                    header_lines += 1
                else:
                    break

        tmp_acc_node = pd.read_csv(assembly_summary,
                                   sep="\t",
                                   header=None,
                                   skiprows=header_lines,
                                   # usecols = 1:assembly_accession, 6:taxid, 8:organism_name, 9:infraspecific_name
                                   usecols=[0, 5, 7, 8],
                                   names=["target", "node", "organism_name", "infraspecific_name"],
                                   index_col="target",
                                   converters={"target": lambda x: x if x in unique_acc else None, "node": str})
        tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()]  # keep only seqids used

        # save count to return
        count_assembly_summary[assembly_summary] = tmp_acc_node.shape[0]
        if not count_assembly_summary[assembly_summary]:
            continue

        # Create specialization
        if level == "assembly":
            # infraspecific_name has a prefix: breed=, cultivar=, ecotype= or strain=
            tmp_acc_node["infraspecific_name"] = tmp_acc_node["infraspecific_name"].replace("^[a-z]+=", "", regex=True).fillna("")

            def build_name(n):
                if n.organism_name.endswith(n.infraspecific_name):
                    return n.organism_name
                else:
                    return n.organism_name + " " + n.infraspecific_name

            # add sufix of the infraspecific_name if not yet contained in the end of organism_name
            tmp_acc_node["specialization_name"] = tmp_acc_node[["organism_name",
                                                                "infraspecific_name"]].apply(lambda n: build_name(n),
                                                                                             axis=1)
            tmp_acc_node["specialization"] = tmp_acc_node.index

        # merge node(taxid) and specialization retrieved based on target(accesion)
        if count_assembly_summary[assembly_summary]:
            info.update(tmp_acc_node)
        del tmp_acc_node

        #if already found all seqids no need to parse all files till the end)
        if sum(count_assembly_summary.values()) == len(unique_acc):
            break

    return count_assembly_summary


def run_eutils(cfg, info, build_output_folder, skip_taxid: bool=False, level: str=""):
    """
    run ganon-get-seq-info.sh script to retrieve taxonomic and assembly info from sequence accessions (ncbi)
    """
    # Write target (index) to a file
    accessions_file = build_output_folder + "accessions.txt"
    rm_files(accessions_file)
    info.to_csv(accessions_file, columns=[], header=False)

    # (-k) always return all entries in the same order
    # (-e) get taxid length
    # (-a) get assembly accession
    # (-m) get assembly name
    # || true to ignore exit status in case some sequences were not retrieved
    run_get_seq_info_cmd = "{0} -i {1} -k {2} {3} || true".format(cfg.path_exec["get_seq_info"],
                                                                  accessions_file,
                                                                  "" if skip_taxid else "-e",
                                                                  "-a -m" if level == "assembly" else "")

    stdout = run(run_get_seq_info_cmd, ret_stdout=True, shell=True, quiet=cfg.quiet)

    # set "na" as NaN with na_values="na"
    if level == "assembly":
        # return target, [taxid,] specialization, specialization_name
        if skip_taxid:
            return pd.read_csv(StringIO(stdout),
                               sep="\t",
                               names=["target", "specialization", "specialization_name"],
                               index_col="target",
                               header=None,
                               dtype=object,
                               na_values="na")
        else:
            return pd.read_csv(StringIO(stdout),
                               sep="\t",
                               names=["target", "length", "node", "specialization", "specialization_name"],
                               index_col="target",
                               header=None,
                               usecols=["target", "node", "specialization", "specialization_name"],
                               dtype=object,
                               na_values="na")
    else:
        # return target, taxid
        return pd.read_csv(StringIO(stdout),
                           sep="\t",
                           names=["target", "length", "node"],
                           index_col="target",
                           header=None,
                           usecols=["target", "node"],
                           dtype=object,
                           na_values="na")
