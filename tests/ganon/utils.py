import shutil
import os
import gzip
import sys
import pandas as pd
from pathlib import Path
from math import floor

from multitax import CustomTx

sys.path.append('src')
from ganon.util import download
from ganon.config import Config
from ganon import ganon

def run_ganon(cfg, prefix):
    """
    capture stderr to log file and run
    """
    with open(prefix + ".log", "w") as sys.stderr:
        return ganon.main(cfg=cfg)


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


def list_files_folder(folder, ext: str = "", recursive: bool = False):
    file_list = []
    if recursive:
        for path in Path(folder).rglob('*' + ext):
            file_list.append(str(path.joinpath()))
    else:
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
    colums = ['file', 'target']
    types = {'file': 'str', 'target': 'str'}
    return pd.read_table(target_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_info(info_file):
    colums = ['file', 'target', 'node',
              'specialization', 'specialization_name']
    types = {'file': 'str', 'target': 'str', 'node': 'str',
             'specialization': 'str', 'specialization_name': 'str'}
    return pd.read_table(info_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_tre(tre_file, output_format: str = "tsv"):
    colums = ['rank', 'target', 'lineage', 'name', 'unique',
              'shared', 'children', 'cumulative', 'cumulative_perc']
    types = {'rank': 'str', 'target': 'str', 'lineage': 'str', 'name': 'str', 'unique': 'uint64',
             'shared': 'uint64', 'children': 'uint64', 'cumulative': 'uint64', 'cumulative_perc': 'float'}
    return pd.read_table(tre_file, sep="," if output_format == "csv" else "\t", header=None, skiprows=0, names=colums, dtype=types)


def parse_tsv(tsv_file):
    return pd.read_table(tsv_file, sep='\t', index_col=0)


def parse_all_one(file):
    colums = ['readid', 'target', 'count']
    types = {'readid': 'str', 'target': 'str', 'count': 'uint64'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def parse_rep(rep_file):
    colums = ['hierarchy', 'target', 'total', 'unique', 'lca', 'rank', 'name']
    types = {'hierarchy': 'str', 'target': 'str', 'total': 'str',
             'unique': 'str', 'lca': 'str', 'rank': 'str', 'name': 'str'}
    return pd.read_table(rep_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)


def build_sanity_check_and_parse(params, skipped_targets: bool = False):
    # Provide sanity checks for outputs (not specific to a test) and returns loaded data
    res = {}

    if not check_files(params["db_prefix"], ["ibf"] if params["taxonomy"] == "skip" else ["ibf", "tax"]):
        if not check_files(params["db_prefix"], ["hibf"] if params["taxonomy"] == "skip" else ["hibf", "tax"]):
            return None

    # target_info file to be read by ganon-build
    # res["target"] fields ['file', 'target']
    res["target"] = parse_target(
        params["db_prefix"] + "_files/build/target_info.tsv")

    if not skipped_targets:
        # Count number of target based on input files
        ntarget = 0
        if params["input_file"]:
            with open(params["input_file"], "r") as fp:
                ntarget = len(fp.readlines())
        else:
            input_files = []
            for i in params["input"]:
                if os.path.isdir(i):
                    input_files.extend(list_files_folder(
                        i, ext=params["input_extension"], recursive=params["input_recursive"]))
                else:
                    input_files.append(i)

            if params["input_target"] == "sequence":
                ntarget = len(list_sequences(input_files))
            else:
                ntarget = len(input_files)

        # Correct number of targets
        if res["target"].shape[0] != ntarget:
            print("Wrong number of targets")
            return None

    # res["info"] fields ['file', 'target', 'node', 'specialization', 'specialization_name']
    if "input_file" in params and params["input_file"]:
        # Parse info file if provided
        res["info"] = parse_info(params["input_file"])
        if res["info"]["target"].isna().all(): # in case of one col only, use filename as target
            res["info"]["target"] = res["info"]["file"].apply(lambda x: Path(x).name)
    else:
        # Generated on build process
        res["info"] = parse_info(params["db_prefix"] + ".info.tsv")

    # res["tax"] = CustomTx from multitax
    if params["taxonomy"] != "skip":
        # Load with multitax, verify consistency
        try:
            res["tax"] = CustomTx(
                files=params["db_prefix"] + ".tax", cols=["node", "parent", "rank", "name"])
        except AssertionError:
            print("Inconsistent .tax")
            return None

        # All targets should be in the tax
        if any(res["target"]["target"].apply(res["tax"].parent) == res["tax"].undefined_node):
            print("Target entries missing on .tax")
            return None

    # Validate specialization and targets
    if params["level"] in ["assembly", "custom"]:  # specializations
        # Check if target_info has the correct value
        if not (res["target"]["target"] == res["info"]["specialization"]).all():
            print("Wrong target")
            return None

        # Check if specialization rank is present in tax
        if params["taxonomy"] != "skip":
            if params["level"] not in res["tax"]._ranks.values():
                print("Rank missing from tax")
                return None
    else:
        # Check if no specialization was generated
        if not res["info"]["specialization"].isna().all() or not res["info"]["specialization_name"].isna().all():
            print("Unrequested specialization")
            return None

        if params["level"]:  # If given and not assembly or custom, should be a taxonomic rank or leaves
            # Check if specialization rank is present in tax (skip for leaves, not a specific rank)
            if params["taxonomy"] != "skip" and params["level"] != "leaves":
                if params["level"] not in res["tax"]._ranks.values():
                    print("Rank missing from tax")
                    return None
            else:
                if not (res["target"]["target"] == res["info"]["node"]).all():
                    print("Wrong target")
                    return None

        else:  # if not given, default to --input-target
            if not (res["target"]["target"] == res["info"]["target"]).all():
                print("Wrong target")
                return None

            # Check if specialization rank is present in tax
            if params["taxonomy"] != "skip":
                if params["input_target"] not in res["tax"]._ranks.values():
                    print("Rank missing from tax")
                    return None

    return res


def classify_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    out_ext = ["rep"]
    if not params["skip_report"]:
        out_ext.append("tre")

    if params["output_one"]:
        out_ext.append("one")

    if params["output_all"]:
        out_ext.append("all")

    if not check_files(params["output_prefix"], out_ext):
        return None

    res = {}
    if params["output_all"]:
        res["all_pd"] = parse_all_one(params["output_prefix"]+".all")
        if res["all_pd"].empty:
            return None

    if params["output_one"]:
        res["one_pd"] = parse_all_one(params["output_prefix"]+".one")
        if res["one_pd"].empty:
            return None

    if "tre" in out_ext:
        res["tre_pd"] = parse_tre(params["output_prefix"]+".tre")
        if res["tre_pd"].empty:
            return None

    res["rep_pd"] = parse_rep(params["output_prefix"]+".rep")
    if res["rep_pd"].empty:
        return None

    return res


def reassign_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    res = {}
    if not params["remove_all"]:
        res["all_pd"] = parse_all_one(params["input_prefix"]+".all")
        if res["all_pd"].empty:
            return None

    if not params["skip_one"]:
        # If no .one file was created
        if not check_files(params["output_prefix"], ["one"]):
            return None
        else:
            res["one_pd"] = parse_all_one(params["output_prefix"]+".one")
            if res["one_pd"].empty:
                return None

    if check_files(params["output_prefix"], ["rep"]):
        res["rep_pd"] = parse_rep(params["output_prefix"] + ".rep")
        if res["rep_pd"].empty:
            return None

    return res


def report_sanity_check_and_parse(params, sum_full_percentage: bool = True):

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
        res["tre_pd"] = parse_tre(
            out_tre, params["output_format"] if "output_format" in params else None)

        # get idx for root (idx_root) and root + unclassified (idx_base)
        # strip white spaces for output_format text
        res["idx_root"] = res["tre_pd"]['rank'].str.strip() == "root"
        if params["report_type"] == "matches":
            res["idx_base"] = res["idx_root"]
        else:
            res["idx_base"] = res["idx_root"] | (
                res["tre_pd"]['rank'] == "unclassified")

        # Check if total is 100%
        if sum_full_percentage and floor(res["tre_pd"][res["idx_base"]]["cumulative_perc"].sum()) != 100:
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

        # check if sum of percentage for each default rank is equal or lower than 100 (floor for rounding errors)
        for rank, val in res["tre_pd"].groupby(by="rank")["cumulative_perc"].sum().apply(floor).items():
            if rank in Config.choices_default_ranks and val > 100:
                print("Inconsistent percentage by rank (>100%)")
                return None

        # Check if counts are consistent
        if ((res["tre_pd"]["unique"] + res["tre_pd"]["shared"] + res["tre_pd"]["children"]) > res["tre_pd"]["cumulative"]).any():
            print("Inconsistent counts")
            return None

        # Check if percentage/abundance is consistent
        # nodes cannot have higher percentage than parents
        target_perc = dict(
            zip(res["tre_pd"]["target"], res["tre_pd"]["cumulative_perc"]))
        for l in res["tre_pd"]["lineage"]:
            # skip empty nodes (e.g 241||412412 -> 241|412412)
            lineage = [n for n in l.split("|") if n]
            # From leaf to root
            for node_idx in list(range(len(lineage)))[::-1]:
                # if node_idx is not latest (0 -> no parent) and if node is reported
                if node_idx and lineage[node_idx] in target_perc:
                    parent = lineage[node_idx-1]
                    # If parent is reported
                    if parent in target_perc and target_perc[lineage[node_idx]] > target_perc[parent]:
                        #print(lineage[node_idx], parent, target_perc[lineage[node_idx]], target_perc[parent])
                        print("Inconsistent percentage among children")

        multi_res[out_tre] = res

    # If only one output, return directly
    if len(multi_res) == 1:
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
    if not (res["out_pd"].values >= 0).all():
        print("Invalid negative value")
        return None

    # specific for percentage output
    if params["output_value"] == "percentage":
        # check if all percentages are on the 0-1 range
        if not (res["out_pd"].values <= 1).all():
            print("Invalid table percentage (>1)")
            return None

        # Check if sum of the lines do not pass 1
        if not (res["out_pd"].sum(axis=1) <= 1).all():
            print("Invalid line value (>1)")
            return None

    return res


def download_bulk_files(download_dir):

    bulk_files = {"https://ftp.ncbi.nlm.nih.gov/":
                  ["pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz",
                      "pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                      "genomes/refseq/assembly_summary_refseq.txt",
                      "genomes/refseq/assembly_summary_refseq_historical.txt",
                      "genomes/genbank/assembly_summary_genbank.txt",
                      "genomes/genbank/assembly_summary_genbank_historical.txt",
                      "genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz"],
                  "https://data.gtdb.ecogenomic.org/releases/latest/":
                      ["ar53_metadata.tsv.gz",
                       "bac120_metadata.tsv.gz"]
                  }

    for url, files in bulk_files.items():
        for file in files:
            if os.path.isfile(download_dir + file):
                print("File found: " + download_dir + file)
            else:
                p = os.path.dirname(file)
                print("Downloading " + url + file +
                      " -> " + download_dir + p + file)
                os.makedirs(download_dir + p, exist_ok=True)
                download([url + file], download_dir + p + "/")


def write_input_file(files,
                     assembly_summary,
                     output_file,
                     cols: list = ["file", "target", "node", "specialization", "specialization_name"],
                     sequence: bool = False):
    """
    Generates a simulated --input-file based on an assembly summary and a list of files
    """

    input_file = {}
    for file in files:
        input_file[Path(file).name] = {"file": file}

    with open(assembly_summary, "r") as asf:
        for line in asf:
            fields = line.split("\t")
            # filename from url
            file = Path(fields[19]).name + "_genomic.fna.gz"

            if file not in input_file:
                continue

            # assembly accession as target
            target = fields[0]
            # taxid as node
            node = fields[5]
            # use biosample as specialization
            spec = fields[2]
            # use name + strain for specialization name
            spec_name = fields[7] + " " + fields[8]

            input_file[file].update({"target": target, "node": node, "specialization": spec, "specialization_name": spec_name})

    with open(output_file, "w") as outf:
        for file, vals in input_file.items():
            if sequence:
                for seq in list_sequences([vals["file"]]):
                    vals["target"] = seq
                    print("\t".join([v for k,v in vals.items() if k in cols]), file=outf)    
            else:
                print("\t".join([v for k,v in vals.items() if k in cols]), file=outf)
