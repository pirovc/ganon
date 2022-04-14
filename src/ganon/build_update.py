import time
import pandas as pd
import re
import os

from ganon.util import validate_input_files
from ganon.util import print_log
from ganon.util import set_tmp_folder
from ganon.util import run
from ganon.util import check_file
from ganon.util import download
from io import StringIO

from multitax import NcbiTx, GtdbTx

def build(cfg):
    return True

def update(cfg):
    return True

def build_custom(cfg):

    # Retrieve and check input files or folders
    input_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet)
    if not input_files:
        print_log("ERROR: No valid input files found")
        return False
    print_log("")

    # Set --input-target if not manually set
    if not cfg.input_target:
        cfg.input_target = "sequence" if len(input_files)==1 else "file"

    # Set working folder 
    tmp_output_folder = cfg.db_prefix + "_tmp/"
    if not set_tmp_folder(tmp_output_folder, cfg.restart): return False

    # Set-up taxonomy
    if cfg.taxonomy:
        tx = time.time()
        print_log("Parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)
        if cfg.taxonomy=="ncbi":
            tax = NcbiTx(files=cfg.taxonomy_files)
        elif cfg.taxonomy=="gtdb":
            tax = GtdbTx(files=cfg.taxonomy_files)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n\n", cfg.quiet)
    else:
        tax = None

    # Set-up input info
    info = load_input(cfg, input_files)

    # Retrieve target info (if not provided as file)
    if not cfg.input_file:
        if cfg.input_target=="sequence":
            get_sequence_info(cfg, info, tmp_output_folder)
        else:
            get_file_info(cfg, info, tmp_output_folder)

    print(info)

    # # If taxonomy is given and level is an taxonomic rank, convert nodes
    # if tax:
    #     # get latest tax.latest()
    #     # validate node on given tax (if came from target_info for example and it's outdated or wrong)
    #     if cfg.level!="leaves" and cfg.level in set(tax._ranks):
    #         info = replace_node_rank()

    # else:
    #     # na on col nodes

    # replaced_spec = validate_specialization(info)
    # if replaced_spec:
    #     print_log(str(replaced_spec) + " invalid specialization entries replaced by target\n", cfg.quiet)

    # filter tax
    # write tax

    # write target-info file

    # run ganon-build

    return True

def update_custom(cfg):
    return True


def load_input(cfg, input_files):
    """
    Load basic target info, either provided as --target-info
    or extracted from file/sequences
    """
    target_info_colums = ["target", "node", "specialization", "file"]
    tx = time.time()
    if cfg.input_file:
        print_log("Parsing --target-info " + cfg.input_file, cfg.quiet)
        info = pd.read_csv(cfg.input_file, sep="\t", header=None, skiprows=0, index_col=None, dtype="str", names=target_info_colums)
    else:
        if cfg.input_target=="sequence":
            print_log("Parsing sequences from --input (" + str(len(input_files)) + " files)", cfg.quiet)
            info = parse_sequence_accession(input_files, target_info_colums)
        else:
            print_log("Parsing --input (" + str(len(input_files)) + " files)", cfg.quiet)
            info = parse_file_accession(input_files, target_info_colums)

    # Drop cols without values
    info.dropna(how="all", inplace=True)
    # Drop duplicated target (seqid or fileid)
    info.drop_duplicates(subset=['target'], inplace=True)
    # set target as index
    info.set_index('target', inplace=True)
    print_log(" - " + str(info.shape[0]) + " unique entries", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return info
        
def validate_specialization(info):
    """
    validate specializatio, if given
    each specialization can have only one parent node
    return number of replaced specializations
    """

    # check for invalid specialization entries
    idx_null_spec = info.specialization.isnull()
    # if all entries are null, no specialization was given
    if all(idx_null_spec):
        return 0
    # get unique tuples node-specialization
    node_spec = info[['node', 'specialization']].drop_duplicates()
    # check for duplicated specialization in the tuples
    idx_multi_parent_spec = info.specialization.isin(node_spec.specialization[node_spec.specialization.duplicated(keep=False)].unique())
    # merge indices for invalid entries
    idx_replace = idx_null_spec | idx_multi_parent_spec
    if idx_replace.any():
        # replace invalid specialization entries with target
        info.loc[idx_replace,"specialization"] = info.index[idx_replace]
        return sum(idx_replace)
    return 0

def parse_sequence_accession(input_files, target_info_colums):
    """
    Look for sequence accession (anything from > to the first space) in all input files
    """
    info = pd.DataFrame(columns=target_info_colums)

    for file in input_files:
        # cat | zcat | gawk -> compability with osx
        run_get = "cat {0} {1} | grep -o '^>[^ ]*' | sed 's/>//'".format(file, "| zcat" if file.endswith(".gz") else "")
        stdout, stderr = run(run_get, shell=True)
        tmp_info = pd.read_csv(StringIO(stdout), header=None, names=['target'])
        tmp_info["file"] = file
        info = pd.concat([info, tmp_info])
    return info


def parse_file_accession(input_files, target_info_colums):
    """
    Look for genbank/refseq assembly accession* pattern in the filename
    if not found, return basename of the file as target

    *https://support.nlm.nih.gov/knowledgebase/article/KA-03451/en-us
    *https://https.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/gca-and-gcf-explained/
    """
    info = pd.DataFrame(columns=target_info_colums)
    assembly_accessions = []
    assembly_accession_pattern = re.compile("GC[A|F]_[0-9]+\.[0-9]+")
    for file in input_files:
        match = assembly_accession_pattern.search(file)
        assembly_accessions.append((match.group() if match else os.path.basename(file), file))
    info[["target","file"]] = assembly_accessions
    return info


def get_file_info(cfg, info, tmp_output_folder):
    assembly_summary_urls = []
    assembly_summary_files = []

    for assembly_summary in cfg.get_file_info:
        # Either a file or prefix for download
        if assembly_summary in cfg.choices_get_file_info:
            assembly_summary_urls.append("ftp://ftp.ncbi.nlm.nih.gov/genomes/" + assembly_summary + "/assembly_summary_" + assembly_summary + ".txt")
        else:
            assembly_summary_files.append(assembly_summary)
        
    if assembly_summary_urls:
        tx = time.time()
        print_log("Downloading assembly_summary files", cfg.quiet)
        assembly_summary_files.extend(download(assembly_summary_urls, tmp_output_folder))
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    tx = time.time()
    print_log("Parsing assembly_summary files", cfg.quiet)
    count_assembly_summary = parse_assembly_summary(info, assembly_summary_files, cfg.level)
    for assembly_summary_file, cnt in count_assembly_summary.items():
        print_log(" - " + str(cnt) + " entries found in the " + assembly_summary_file.split("/")[-1] + " file", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)


def get_sequence_info(cfg, info, tmp_output_folder):
    # Max. sequences to use eutils in auto mode
    max_seqs_eutils = 50000

    if not cfg.get_sequence_info:
        # auto
        if info.shape[0] > max_seqs_eutils: 
            mode = ["nucl_gb", "nucl_wgs"]
        else:
            mode = ["eutils"]
    elif "eutils" in cfg.get_sequence_info:
        mode = ["eutils"]
    else:
        mode = cfg.get_sequence_info # custom accession2taxid prefix or files

    if mode[0]=="eutils":
        tx = time.time()
        print_log("Retrieving sequence information from NCBI e-utils", cfg.quiet)
        run_eutils(cfg, info, tmp_output_folder)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        acc2txid_urls = []
        acc2txid_files = []
        for acc2txid in mode:
            # Either a file or prefix for download
            if acc2txid in cfg.choices_get_sequence_info:
                acc2txid_urls.append("ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/" + acc2txid + ".accession2taxid.gz")
            else:
                acc2txid_files.append(acc2txid)                

        if acc2txid_urls:
            tx = time.time()
            print_log("Downloading accession2taxid files", cfg.quiet)
            acc2txid_files.extend(download(acc2txid_urls, tmp_output_folder))
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        tx = time.time()
        print_log("Parsing accession2taxid files", cfg.quiet)
        count_acc2txid = parse_acc2txid(info, acc2txid_files)
        for acc2txid_file, cnt in count_acc2txid.items():
            print_log(" - " + str(cnt) + " entries found in the " + acc2txid_file.split("/")[-1] + " file", cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        # If assembly name or accession are requested, get with e-utils
        if cfg.level in ["assembly-name", "assembly-acc"]:
            print_log("Retrieving assembly information from NCBI e-utils (--level " + cfg.level + ")", cfg.quiet)
            run_eutils(cfg, info, tmp_output_folder, skip_taxid=True)
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
                                   converters={"target": lambda x: x if x in unique_acc else None},
                                   dtype={"node": "str"})
        tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()] #  keep only seqids used
        tmp_acc_node = tmp_acc_node[tmp_acc_node["node"] != "0"] #  filter out taxid==0

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
                                   skiprows=2,
                                   usecols=[0, 5, 7, 8],  # 1:assembly_accession, 6:taxid, 8:organism_name, 9:infraspecific_name
                                   names=["target", "node", "organism_name", "infraspecific_name"],
                                   index_col="target",
                                   converters={"target": lambda x: x if x in unique_acc else None},
                                   dtype={"node": "str"})
        tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()] #  keep only seqids used
        tmp_acc_node = tmp_acc_node[tmp_acc_node["node"] != "0"] #  filter out taxid==0

        # Create specialization if requested for --level
        if level == "assembly-name":
            tmp_acc_node["specialization"] = tmp_acc_node["organism_name"] + " " + tmp_acc_node["infraspecific_name"]
        elif level == "assembly-acc":
            tmp_acc_node["specialization"] = tmp_acc_node.index

        # save count to return
        count_assembly_summary[assembly_summary] = tmp_acc_node.shape[0]
        # merge node(taxid) and specialization retrieved based on target(accesion)
        if count_assembly_summary[assembly_summary]:
            info.update(tmp_acc_node)
        del tmp_acc_node

        #if already found all seqids no need to parse all files till the end)
        if sum(count_assembly_summary.values()) == len(unique_acc):
            break

    return count_assembly_summary


def run_eutils(cfg, info, tmp_output_folder, skip_taxid: bool = False):
    """
    run ganon-get-seq-info.sh script to retrieve taxonomic and assembly info from sequence accessions (ncbi)
    """
    # Write target (index) to a file
    accessions_file = tmp_output_folder + "accessions.txt"
    info.to_csv(accessions_file, columns=[], header=False)

    get_assembly = ""
    if cfg.level == "assembly-acc":
        get_assembly = "-a"
    elif cfg.level == "assembly-name":
        get_assembly = "-m"

    # (-k) always return all entries in the same order 
    # (-r) report sequence accession if assembly not found (-r)
    run_get_seq_info_cmd = "{0} -i {1} -k -r {2} {3}".format(cfg.path_exec["get_seq_info"],
                                                            accessions_file,
                                                            "" if skip_taxid else "-e",
                                                            get_assembly)
    stdout, stderr = run(run_get_seq_info_cmd, print_stderr=True if not cfg.quiet else False, exit_on_error=False)

    # set "na" as NaN with na_values="na"
    if get_assembly:
        if skip_taxid:
            info.update(pd.read_csv(StringIO(stdout), sep="\t", names=["target", "specialization"], index_col="target", header=None, skiprows=0, dtype="str", na_values="na"))
        else:
            info.update(pd.read_csv(StringIO(stdout), sep="\t", names=["target","length","node","specialization"], index_col="target", header=None, skiprows=0, usecols=["target", "node", "specialization"], dtype="str", na_values="na"))
    else:
        info.update(pd.read_csv(StringIO(stdout), sep="\t", names=["target","length","node"], index_col="target", header=None, skiprows=0, usecols=["target", "node"], dtype="str", na_values="na"))

