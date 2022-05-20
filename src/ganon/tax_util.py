import time
import pandas as pd
import re
import os

from ganon.util import download, run, print_log, rm_files
from io import StringIO


def parse_sequence_accession(input_files, info):
    """
    Look for sequence accession (anything from > to the first space) in all input files
    """
    for file in input_files:
        # cat | zcat | gawk -> compability with osx
        run_cat = "cat {0} {1} | grep -o '^>[^ ]*' | sed 's/>//'".format(file, "| zcat" if file.endswith(".gz") else "")
        stdout = run(run_cat, shell=True, ret_stdout=True)
        tmp_info = pd.read_csv(StringIO(stdout), header=None, names=['target'])
        tmp_info["file"] = file
        info = pd.concat([info, tmp_info])
    return info


def parse_file_accession(input_files, info):
    """
    Look for genbank/refseq assembly accession* pattern in the filename
    if not found, return basename of the file as target

    *https://support.nlm.nih.gov/knowledgebase/article/KA-03451/en-us
    *https://https.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/gca-and-gcf-explained/
    """
    assembly_accessions = []
    assembly_accession_pattern = re.compile("GC[A|F]_[0-9]+\.[0-9]+")
    for file in input_files:
        match = assembly_accession_pattern.search(file)
        assembly_accessions.append((match.group() if match else os.path.basename(file), file))
    info[["target", "file"]] = assembly_accessions
    return info


def get_file_info(cfg, info, tax, build_output_folder):
    if cfg.taxonomy == "ncbi" or (cfg.taxonomy == "skip" and cfg.level == "assembly"):
        assembly_summary_urls = []
        assembly_summary_files = []

        for assembly_summary in cfg.ncbi_file_info:
            # Either a file or prefix for download
            if assembly_summary in cfg.choices_ncbi_file_info:
                # split if _historical
                assembly_summary_urls.append(cfg.ncbi_ftp +
                                             "/genomes/" + assembly_summary.split("_")[0] + 
                                             "/assembly_summary_" + assembly_summary + ".txt")
            else:
                assembly_summary_files.append(assembly_summary)

        if assembly_summary_urls:
            tx = time.time()
            print_log("Downloading assembly_summary files", cfg.quiet)
            assembly_summary_files.extend(download(assembly_summary_urls, build_output_folder))
            # Skip headers (first 2 lines) so pandas.read_csv read it properly
            for f in assembly_summary_files:
                with open(f, 'r+') as fp:
                    lines = fp.readlines()
                    fp.seek(0)
                    fp.truncate()
                    fp.writelines(lines[2:])
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
                        acc2txid_urls.append(cfg.ncbi_ftp + "/pub/taxonomy/accession2taxid/" +
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
    for acc2txid in acc2txid_files:
        tmp_acc_node = pd.read_csv(acc2txid,
                                   sep="\t",
                                   header=None,
                                   skiprows=1,
                                   usecols=[1, 2],
                                   names=["target", "node"],
                                   index_col="target",
                                   converters={"target": lambda x: x if x in unique_acc else None, "node": str})
        tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()]  # keep only seqids used
        tmp_acc_node = tmp_acc_node[tmp_acc_node["node"] != "0"]  # filter out taxid==0

        # save count to return
        count_acc2txid[acc2txid] = tmp_acc_node.shape[0]
        # merge node(taxid) retrieved based on target(accesion)
        if count_acc2txid[acc2txid]:
            info.update(tmp_acc_node)
        del tmp_acc_node
        #if already found all seqids no need to parse all files till the end)
        if sum(count_acc2txid.values()) == len(unique_acc):
            break

    return count_acc2txid


def parse_assembly_summary(info, assembly_summary_files, level):
    count_assembly_summary = {}
    unique_acc = set(info.index)

    for assembly_summary in assembly_summary_files:
        tmp_acc_node = pd.read_csv(assembly_summary,
                                   sep="\t",
                                   header=None,
                                   # skiprows=2, do not skiprows (no header from genome_updater or skipped before)
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
    run_get_seq_info_cmd = "{0} -i {1} -k {2} {3}".format(cfg.path_exec["get_seq_info"],
                                                          accessions_file,
                                                          "" if skip_taxid else "-e",
                                                          "-a -m" if level == "assembly" else "")
    stdout = run(run_get_seq_info_cmd, ret_stdout=True, quiet=cfg.quiet)

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
