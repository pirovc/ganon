import shutil, os
import pandas as pd
from pathlib import Path
from ganon.bins import Bins
from ganon.gnn import Gnn

def setup_dir(results_dir):
    shutil.rmtree(results_dir, ignore_errors=True)
    os.makedirs(results_dir)

def check_files(prefix, extensions):
    for ext in extensions:
        f=prefix+"."+ext if ext else prefix
        if not Path(f).is_file():
            print("File (" + f +") was not created")
            return False
        elif Path(f).stat().st_size==0:
            print("File (" + f +") is empty")
            return False
    return True
    
def parse_bins(bins):
    #columns=['seqid', 'seqstart', 'seqend', 'length', 'taxid', 'binid', 'specialization']
    #types={'seqid': 'str', 'seqstart': 'uint64', 'seqend': 'uint64', 'length': 'uint64', 'taxid': 'str', 'binid': 'uint64', 'specialization': 'str'}
    bins.bins['binid'] = pd.to_numeric(bins.bins['binid'])
    bins.bins['length'] = pd.to_numeric(bins.bins['length'])
    return bins.bins

def parse_seq_info(seq_info_file):
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    return pd.read_table(seq_info_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_map(map_file):
    colums=['target', 'binid']
    types={'target': 'str', 'binid': 'uint64'}
    return pd.read_table(map_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_tax(tax_file):
    colums=['taxid', 'parent', 'rank', 'name']
    types={'taxid': 'str', 'parent': 'str', 'rank': 'str', 'name': 'str'}
    return pd.read_table(tax_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_tre(tre_file):
    colums=['rank', 'target', 'lineage', 'name', 'unique', 'total', 'cumulative', 'cumulative_perc']
    types={'rank': 'str', 'target': 'str', 'lineage': 'str', 'name': 'str', 'unique': 'str', 'total': 'str', 'cumulative': 'str', 'cumulative_perc': 'str'}
    return pd.read_table(tre_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_tsv(tsv_file):
    return pd.read_table(tsv_file, sep='\t', index_col=0)

def parse_all_lca(file):
    colums=['readid', 'target', 'count']
    types={'readid': 'str', 'target': 'str', 'count': 'uint64'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_rep(rep_file):
    colums=['hierarchy', 'target', 'total', 'unique', 'lca', 'rank', 'name']
    types={'hierarchy': 'str', 'target': 'str', 'total': 'str', 'unique': 'str', 'lca': 'str', 'rank': 'str', 'name': 'str'}
    return pd.read_table(rep_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def build_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["db_prefix"], ["ibf", "map", "tax", "gnn"]):
        return None

    res = {}
    # Parse in and out files
    if "seq_info_file" in params and params["seq_info_file"]:
        res["seq_info"] = parse_seq_info(params["seq_info_file"])
    else:
        res["seq_info"] = parse_seq_info(params["db_prefix"]+".seqinfo.txt")
         
    res["gnn"] = Gnn(file=params["db_prefix"]+".gnn")
    res["tax_pd"] = parse_tax(params["db_prefix"]+".tax")
    res["map_pd"] = parse_map(params["db_prefix"]+".map")
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

def report_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_prefix"], ["tre"]):
        return None

    res = {}
    # Sequence information from database to be updated
    res["tre_pd"] =  parse_tre(params["output_prefix"] + ".tre")
    return res

def table_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_file"], [""]):
        return None

    res = {}
    # Sequence information from database to be updated
    res["tsv_pd"] = parse_tsv(params["output_file"])
    return res

def update_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_db_prefix"], ["ibf", "map", "tax", "gnn"]):
        return None

    res = {}
    # Sequence information from database to be updated
    if not params["update_complete"]:
        res["seq_info"] = parse_seq_info(params["db_prefix"]+".seqinfo.txt")
    else: # Do not load it in case of update_complete, where all sequences must be provided
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
