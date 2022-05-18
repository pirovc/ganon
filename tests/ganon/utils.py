import shutil
import os
import gzip
import pandas as pd
from pathlib import Path
from math import floor

from multitax import CustomTx


def setup_dir(results_dir):
    shutil.rmtree(results_dir, ignore_errors=True)
    os.makedirs(results_dir)


def check_files(prefix, extensions):
    for ext in extensions:
        f = prefix + "." + ext if ext else prefix
        if not Path(f).is_file():
            print("File (" + f + ") was not created")
            return False
        elif Path(f).stat().st_size == 0:
            print("File (" + f + ") is empty")
            return False
    return True


def list_files_folder(folder, ext: str=None):
    file_list = []
    for file in os.listdir(folder):
        if ext is None or file.endswith(ext):
            file_list.append(folder + file)
    return file_list


def list_sequences(files):
    sequence_list = []
    for file in files:
        if file.endswith(".gz"):
            f = gzip.open(file, "rt")
        else:
            f = open(file, "r")
        for line in f:
            if line[0] == ">":
                sequence_list.append(line.split(" ")[0][1:])
        f.close()
    return sequence_list


def parse_target(target_file):
    colums = ['file', 'target', 'sequence']
    types = {'file': 'str', 'target': 'str', 'sequence': 'str'}
    return pd.read_table(target_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_info(info_file):
    colums = ['file', 'target', 'node', 'specialization', 'specialization_name']
    types = {'file': 'str', 'target': 'str', 'node': 'str', 'specialization': 'str', 'specialization_name': 'str'}
    return pd.read_table(info_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_tre(tre_file):
    colums = ['rank', 'target', 'lineage', 'name', 'unique', 'shared', 'children', 'cumulative', 'cumulative_perc']
    types = {'rank': 'str', 'target': 'str', 'lineage': 'str', 'name': 'str', 'unique': 'uint64', 'shared': 'uint64', 'children': 'uint64', 'cumulative': 'uint64', 'cumulative_perc': 'float'}
    return pd.read_table(tre_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_tsv(tsv_file):
    return pd.read_table(tsv_file, sep='\t', index_col=0)


def parse_all_lca(file):
    colums = ['readid', 'target', 'count']
    types = {'readid': 'str', 'target': 'str', 'count': 'uint64'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_rep(rep_file):
    colums = ['hierarchy', 'target', 'total', 'unique', 'lca', 'rank', 'name']
    types = {'hierarchy': 'str', 'target': 'str', 'total': 'str', 'unique': 'str', 'lca': 'str', 'rank': 'str', 'name': 'str'}
    return pd.read_table(rep_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def build_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and returns loaded data
    res = {}

    if not check_files(params["db_prefix"], ["ibf"] if params["taxonomy"] == "skip" else ["ibf", "tax"]):
        return None

    # target_info file to be read by ganon-build
    # res["target"] fields ['file', 'target', 'sequence']
    res["target"] = parse_target(params["db_prefix"] + "_files/build/target_info.tsv")

    # res["info"] fields ['file', 'target', 'node', 'specialization', 'specialization_name']
    if "input_file" in params and params["input_file"]:
        # Parse info file if provided
        res["info"] = parse_info(params["input_file"])
    else:
        # Generated on build process
        res["info"] = parse_info(params["db_prefix"] + ".info.tsv")

    # res["tax"] = CustomTx from multitax
    if params["taxonomy"] != "skip":
        # Load with multitax, verify consistency
        try:
            res["tax"] = CustomTx(files=params["db_prefix"] + ".tax", cols=["node", "parent", "rank", "name"])
        except AssertionError:
            print("Inconsistent .tax")
            return None

        # All targets should be in the tax
        if any(res["target"]["target"].apply(res["tax"].parent) == res["tax"].undefined_node):
            print("Target entries missing on .tax")
            return None

    return res


def classify_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    out_ext = ["lca","rep","tre"]
    if params["output_all"]: out_ext.append("all")
    if not check_files(params["output_prefix"], out_ext):
        return None

    res = {}
    if params["output_all"]:
        res["all_pd"] = parse_all_lca(params["output_prefix"]+".all")
        # Sequence information from database to be updated
    res["tre_pd"] = parse_tre(params["output_prefix"]+".tre")
    res["rep_pd"] = parse_rep(params["output_prefix"]+".rep")
    res["lca_pd"] = parse_all_lca(params["output_prefix"]+".lca")
    return res

def report_sanity_check_and_parse(params, sum_full_percentage: bool=True):

    # get all output files referent to the run by name
    output_files = []
    directory = os.path.dirname(params["output_prefix"])
    for file in os.listdir(directory):
        if file.startswith(os.path.splitext(os.path.basename(params["output_prefix"]))[0]) and file.endswith(".tre"):
            output_files.append(os.path.join(directory, file))

    multi_res = {}
    for out_tre in output_files:
        # Provide sanity checks for outputs (not specific to a test) and return loaded data
        if not check_files(out_tre, [""]):
            return None

        res = {}
        # Sequence information from database to be updated
        res["tre_pd"] = parse_tre(out_tre)

        # get idx for root (idx_root) and root + unclassified (idx_base)
        res["idx_root"] = res["tre_pd"]['rank'] == "root"
        if params["report_type"] == "reads":
            res["idx_base"] = res["idx_root"] | (res["tre_pd"]['rank'] == "unclassified")
        else:
            res["idx_base"] = res["idx_root"]

        # Check if total is 100%
        if sum_full_percentage and floor(res["tre_pd"][res["idx_base"]]["cumulative_perc"].sum())!=100:
            print("Inconsistent total percentage")
            return None

        # check if sum of all unique+shared is lower or equal to root and unclassified
        if (res["tre_pd"][~res["idx_base"]]["unique"].sum() + res["tre_pd"][~res["idx_base"]]["shared"].sum()) > res["tre_pd"][res["idx_base"]]["cumulative"].sum():
            print("Inconsistent total counts")
            return None

        # Check if any cumulative_perc is higher than 100
        if (res["tre_pd"]["cumulative_perc"] > 100).any():
            print("Inconsistent percentage (>100%)")
            return None

        # check if sum of percentage for each rank (excluding "no rank") is equal or lower than 100 (floor for rounding errors)
        if(res["tre_pd"][res["tre_pd"]["rank"]!="no rank"].groupby(by="rank")["cumulative_perc"].sum().apply(floor)>100).any():
            print("Inconsistent percentage by rank (>100%)")
            return None

        # Check if counts are consistent
        if ((res["tre_pd"]["unique"] + res["tre_pd"]["shared"] + res["tre_pd"]["children"]) > res["tre_pd"]["cumulative"]).any():
            print("Inconsistent counts")
            return None 

        multi_res[out_tre] = res

    # If only one output, return directly
    if len(multi_res)==1:
        return res
    else:
        return multi_res

def table_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_file"], [""]):
        return None

    res = {}
    # Sequence information from database to be updated
    res["out_pd"] = parse_tsv(params["output_file"])

    # check if all values are positive
    if not (res["out_pd"].values>=0 ).all():
        print("Invalid negative value")
        return None

    # specific for percentage output
    if params["output_value"]=="percentage":
        # check if all percentages are on the 0-1 range
        if not (res["out_pd"].values<=1).all():
            print("Invalid table percentage (>1)")
            return None

        # Check if sum of the lines do not pass 1
        if not (res["out_pd"].sum(axis=1)<=1).all():
            print("Invalid line value (>1)")
            return None


    return res

def update_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_db_prefix"], ["ibf", "map", "tax", "gnn"]):
        return None

    res = {}
    # Sequence information from database to be updated
    if not params["update_complete"]:
        res["seq_info"] = parse_info(params["db_prefix"] + ".info.tsv")
    else: 
        res["seq_info"] = pd.DataFrame()

    # Parse in and out files
    if "seq_info_file" in params and params["seq_info_file"]:
        res["seq_info"] = res["seq_info"].append(parse_seq_info(params["seq_info_file"]), ignore_index=True)
    else:
        res["seq_info"] = res["seq_info"].append(parse_seq_info(params["output_db_prefix"]+".seqinfo.txt"), ignore_index=True)

    res["gnn"] = Gnn(file=params["output_db_prefix"]+".gnn")
    res["tax_pd"] = parse_tax(params["output_db_prefix"]+".tax")
    res["map_pd"] = parse_map(params["output_db_prefix"]+".map")
    res["bins_pd"] = parse_bins(Bins(taxsbp_ret=res["gnn"].bins))

    # Check number of bins
    if res["map_pd"].binid.unique().size != res["gnn"].number_of_bins:
        print("Number of bins do not match between .gnn and .map")
        return None

    # Check if all input accession made it to the bins
    if not res["seq_info"]["seqid"].isin(res["bins_pd"]["seqid"]).all():
        print("Missing sequence accessions on bins")
        return None

    # Check if all taxids/assembly on .map appear on .tax
    if res["tax_pd"]["taxid"].isin(res["map_pd"]["target"].drop_duplicates()).all():
        print("Inconsistent entries between taxonomy (.tax) and bin map (.map)")
        return None

    return res
