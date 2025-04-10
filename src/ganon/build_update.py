import time
import pandas as pd
import os
import shutil
import pickle
import math
from pathlib import Path

from ganon.config import Config
from ganon.util import check_file
from ganon.util import check_folder
from ganon.util import print_log
from ganon.util import run
from ganon.util import validate_input_files
from ganon.util import rm_files
from ganon.util import save_state
from ganon.util import set_output_folder
from ganon.util import load_state
from ganon.tax_util import get_file_info
from ganon.tax_util import get_sequence_info
from ganon.tax_util import parse_sequence_accession
from ganon.tax_util import parse_file_accession
from ganon.tax_util import get_genome_size

from multitax import NcbiTx, GtdbTx


def build(cfg):
    
    # Set paths
    if not cfg.set_paths():
        return False

    files_output_folder = set_output_folder(cfg.db_prefix)
    if cfg.restart:
        restart_build(files_output_folder)

    assembly_summary = files_output_folder + "assembly_summary.txt"

    # Skip if already finished download from previous run
    if load_state("build_download", files_output_folder) and check_file(assembly_summary):
        print_log("Download finished - skipping", cfg.quiet)
    else:
        # If assembly_summary.txt was written and some files were already downloaded, try to fix
        resume_download = False
        if check_file(assembly_summary):
            if check_folder(files_output_folder + get_gu_current_version(assembly_summary) + "/files/"):
                print_log("Incomplete files detected, resuming download\n", cfg.quiet)
                resume_download = True

        tx = time.time()
        print_log("Downloading files from " + ",".join(cfg.source) + " [" + ",".join(cfg.organism_group if cfg.organism_group else cfg.taxid) + "]", cfg.quiet)
        run_genome_updater_cmd = " ".join([cfg.path_exec['genome_updater'],
                                           "-d '" + ",".join(cfg.source) + "'",
                                           "-g '" + ",".join(cfg.organism_group) + "'" if cfg.organism_group else "",
                                           "-T '" + ",".join(cfg.taxid) + "'" if cfg.taxid else "",
                                           "-A " + str(cfg.top) if cfg.top else "",
                                           "-l 'complete genome'" if cfg.complete_genomes else "",
                                           "-c 'reference genome'" if cfg.reference_genomes else "",
                                           "-f 'genomic.fna.gz'",
                                           "-t " + str(cfg.threads),
                                           "-o " + files_output_folder,
                                           "-M " + cfg.taxonomy if cfg.taxonomy=="gtdb" else "",
                                           "-m",
                                           "-N",
                                           "-i" if resume_download else "",
                                           "-s" if cfg.quiet else "",
                                           "-w" if not cfg.verbose else "",
                                           cfg.genome_updater if cfg.genome_updater else ""])
        run(run_genome_updater_cmd, quiet=cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
        save_state("build_download", files_output_folder)

    # get current version from assembly_summary
    input_folder = files_output_folder + get_gu_current_version(assembly_summary) + "/files/"

    build_custom_params = {"input": [input_folder],
                           "input_extension": "fna.gz",
                           "input_recursive": True,
                           "input_target": "file",
                           "ncbi_file_info": [assembly_summary]}

    build_default_params = {"db_prefix": cfg.db_prefix,
                            "level": cfg.level,
                            "taxonomy": cfg.taxonomy,
                            "taxonomy_files": cfg.taxonomy_files,
                            "genome_size_files": cfg.genome_size_files,
                            "threads": cfg.threads,
                            "max_fp": cfg.max_fp,
                            "filter_size": cfg.filter_size,
                            "kmer_size": cfg.kmer_size,
                            "window_size": cfg.window_size,
                            "hash_functions": cfg.hash_functions,
                            "mode": cfg.mode,
                            "min_length": cfg.min_length,
                            "verbose": cfg.verbose,
                            "quiet": cfg.quiet,
                            "ganon_path": cfg.ganon_path,
                            "raptor_path": cfg.raptor_path,
                            "n_refs": cfg.n_refs,
                            "n_batches": cfg.n_batches,
                            "filter_type": cfg.filter_type,
                            "write_info_file": cfg.write_info_file,
                            "keep_files": cfg.keep_files}

    build_custom_params.update(build_default_params)

    build_custom_config = Config("build-custom", **build_custom_params)
    save_config(build_custom_config, files_output_folder + "config.pkl")

    ret_build = build_custom(cfg=build_custom_config,
                             which_call="build")

    if ret_build:
        print_log("", cfg.quiet)
        print_log(files_output_folder + " contains reference sequences and configuration files.", cfg.quiet)
        print_log("Keep this folder if you want to update your database later. Otherwise it can be deleted.", cfg.quiet)
        print_log("", cfg.quiet)

    return ret_build


def update(cfg):

    # Set paths
    if not cfg.set_paths():
        return False

    files_output_folder = set_output_folder(cfg.db_prefix)
    if cfg.restart:
        restart_update(files_output_folder)

    tx = time.time()
    # Skip if already finished download from previous run
    if load_state("update_download", files_output_folder):
        print_log("Download finished - skipping", cfg.quiet)
    else:
        print_log("Downloading updated files", cfg.quiet)
        run_genome_updater_cmd = " ".join([cfg.path_exec['genome_updater'],
                                           "-o " + files_output_folder,
                                           "-m",
                                           "-N",
                                           "-s" if cfg.quiet else "",
                                           "-w" if not cfg.verbose else ""])
        run(run_genome_updater_cmd, quiet=cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
        save_state("update_download", files_output_folder)

    # get current version from assembly_summary
    assembly_summary = files_output_folder + "assembly_summary.txt"
    input_folder = files_output_folder + get_gu_current_version(assembly_summary) + "/files/"

    build_custom_params = {"input": [input_folder],
                           "input_extension": "fna.gz",
                           "input_recursive": True,
                           "input_target": "file",
                           "ncbi_file_info": [assembly_summary]}

    build_default_params = {"db_prefix": cfg.output_db_prefix if cfg.output_db_prefix else cfg.db_prefix,
                            "threads": cfg.threads,
                            "verbose": cfg.verbose,
                            "quiet": cfg.quiet,
                            "ganon_path": cfg.ganon_path,
                            "raptor_path": cfg.raptor_path,
                            "n_refs": cfg.n_refs,
                            "n_batches": cfg.n_batches,
                            "write_info_file": cfg.write_info_file,
                            "keep_files": cfg.keep_files}
    build_custom_params.update(build_default_params)

    loaded_params = load_config(files_output_folder + "config.pkl")
    build_custom_params["level"] = loaded_params["level"]
    build_custom_params["taxonomy"] = loaded_params["taxonomy"]
    build_custom_params["max_fp"] = loaded_params["max_fp"]
    build_custom_params["filter_size"] = loaded_params["filter_size"]
    build_custom_params["kmer_size"] = loaded_params["kmer_size"]
    build_custom_params["window_size"] = loaded_params["window_size"]
    build_custom_params["hash_functions"] = loaded_params["hash_functions"]
    build_custom_params["mode"] = loaded_params["mode"] if "mode" in loaded_params else "avg"  # mode introduce in v1.4.0
    build_custom_params["min_length"] = loaded_params["min_length"] if "min_length" in loaded_params else 0  # mode introduce in v1.6.0    
    # filter_type introduced in v2.0.0, before was --hibf
    if "filter_type" in loaded_params:
        ft = loaded_params["filter_type"]  # current definition
    elif "hibf" in loaded_params:
        ft = "hibf"  # --hibf
    else:
        ft = "ibf"  # default < 2.0.0 was ibf
    build_custom_params["filter_type"] = ft

    build_custom_config = Config("build-custom", **build_custom_params)

    ret_build = build_custom(cfg=build_custom_config,
                             which_call="update")

    if ret_build:
        new_files_output_folder = None
        if cfg.output_db_prefix:
            new_files_output_folder = set_output_folder(cfg.output_db_prefix)
            
        # Change input config to new folder    
        if new_files_output_folder:
            build_custom_config.input = [new_files_output_folder + get_gu_current_version(assembly_summary) + "/files/"]
            build_custom_config.ncbi_file_info = [new_files_output_folder + "assembly_summary.txt"]

        # Save config again (change on db_prefix, input folders)
        save_config(build_custom_config, files_output_folder + "config.pkl")

        # Remove save states from finished update (from base folder)
        clear_states("update", files_output_folder)

        if new_files_output_folder:
            if os.path.isfile(new_files_output_folder + "build/target_info.tsv"):
                # Move target_info.tsv created in the new folder to the old and remove empty folder
                shutil.move(new_files_output_folder + "build/target_info.tsv",
                            files_output_folder + "build/target_info.tsv")
                shutil.rmtree(set_output_folder(cfg.output_db_prefix) + "build/")
            # Move files folder to new output_db_prefix
            os.rename(set_output_folder(cfg.db_prefix), new_files_output_folder)

    return ret_build


def build_custom(cfg, which_call: str="build_custom"):
    
    # Set paths
    if not cfg.set_paths():
        return False

    files_output_folder = set_output_folder(cfg.db_prefix)      # DB_PREFIX_files/
    build_output_folder = files_output_folder + "build/"        # DB_PREFIX_files/build/
    target_info_file = build_output_folder + "target_info.tsv"  # DB_PREFIX_files/build/target_info.tsv

    # calling build_custom internally, already checked folders
    if which_call == "build_custom" and cfg.restart:
        restart_build(files_output_folder)

    # Skip if already finished target_info from previous run
    if load_state(which_call + "_parse", files_output_folder):
        print_log("Parse finished - skipping", cfg.quiet)
    else:
        tax = None
        input_files = []

        # Create tmp build folder if not yet existing
        shutil.rmtree(build_output_folder, ignore_errors=True)
        os.makedirs(build_output_folder, exist_ok=True)

        # Retrieve and check input files or folders
        if cfg.input:
            tx = time.time()
            print_log("Parsing --input", cfg.quiet)
            input_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet, input_recursive=cfg.input_recursive)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
            if not input_files:
                print_log("ERROR: No valid input files found", cfg.quiet)
                return False

        # Set-up taxonomy
        if cfg.taxonomy != "skip":
            tax = load_taxonomy(cfg, build_output_folder)

        # Set-up input info
        info = load_input(cfg, input_files, build_output_folder)
        # Define user bins for writing taxonomy and target info file
        user_bins_col = "target"  # Default as target
        if cfg.level in cfg.choices_level:
            user_bins_col = "specialization"  # if specialization was requested
        elif cfg.level and cfg.level not in cfg.choices_input_target:  # if any other level is provided (leaves, species, ...) and not at sequence of file level
            user_bins_col = "node"

        if info.empty:
            print_log("ERROR: Unable to parse input files", cfg.quiet)
            return False

        # Retrieve target info if taxonomy or specialization is required (and if file is not provided)
        if (tax or cfg.level == "assembly") and not cfg.input_file:
            if cfg.input_target == "sequence":
                get_sequence_info(cfg, info, tax, build_output_folder)
            else:
                get_file_info(cfg, info, tax, build_output_folder)

        # Validate taxonomic node only if taxonomy is provided
        if tax:
            validate_taxonomy(info, tax, cfg)
            if info.empty:
                print_log("ERROR: Unable to match taxonomy to targets", cfg.quiet)
                return False

        # Validate specialization for assembly and custom (required to be after taxonomy)
        if cfg.level in cfg.choices_level:
            validate_specialization(info, cfg.quiet)
            if info.empty:
                print_log("ERROR: Unable to match specialization to targets", cfg.quiet)
                return False

        # Filter and write taxonomy
        if tax:
            unique_nodes = info["node"].unique()

            # Check if targets/specializations are not overlapping with taxids
            if (user_bins_col=="target" and info.index.isin(unique_nodes).any()) or \
               (user_bins_col=="specialization" and info["specialization"].isin(unique_nodes).any()):
                print_log("ERROR: " + user_bins_col + " overlaps with taxonomic identifiers", cfg.quiet)
                return False

            # Get estimates of genome sizes
            genome_sizes = get_genome_size(cfg, unique_nodes, tax, build_output_folder)
            
            # filter only used tax. nodes
            tax.filter(unique_nodes) 

            # write tax with added nodes and genome sizes
            write_tax(cfg.db_prefix + ".tax", info, tax, genome_sizes, user_bins_col, cfg.level, cfg.input_target)

        # If requested, save a copy of the info file to re-run build quicker
        if cfg.write_info_file:
            write_info_file(info, cfg.db_prefix + ".info.tsv")

        # Write aux file for ganon-build
        write_target_info(info, user_bins_col, target_info_file)
        save_state(which_call + "_parse", files_output_folder)

    # Skip if already finished run
    if load_state(which_call + "_run", files_output_folder):
        print_log("Build finished - skipping", cfg.quiet)
    else:

        if cfg.filter_type == "hibf":
            print_log("Building index (raptor)", cfg.quiet)

            tx = time.time()
            print_log("Preparing files", cfg.quiet)
            # rewrite target_info_file with one line for each target
            # symbolic link is created for the first time with the target name
            # this is a workaround to name targets with raptor
            target_files = {}
            with open(target_info_file, 'r') as tif:
                for line in tif:
                    col = line.rstrip().split("\t")
                    if col[1] not in target_files:
                        target_files[col[1]] = []
                    target_files[col[1]].append(col[0])
            n_bins = len(target_files)

            raptor_input_file = build_output_folder + "hibf.txt"
 
            # Create symbolic link with target to raptor
            with open(raptor_input_file, "w") as filehibf:
                for target, files in target_files.items():
                    # raptor v3.0.0 "eats" the . (e.g. GCF_013391805.1 -> GCF_013391805)
                    # raptor v3.0.1 fixes it
                    # raptor v3.0.X "eats" the space (e.g. s__Pectobacterium carotovorum -> s__Pectobacterium)
                    # since input file is a space separated, substitute by placeholder "---", treated at runtime ganon-classify
                    new_target = target.replace(" ", "---")
                    # Select first file
                    first_file = os.path.abspath(files[0])
                    # Get extension(s)
                    suffixes = Path(first_file).suffixes
                    # If last one is gz, get real suffix
                    exts = "".join(suffixes[-2:]) if suffixes[-1]==".gz" else suffixes[-1]
                    target_file = build_output_folder + new_target + exts
                    # Attempt to remove if exists +
                    # Create symbolic link with correct name for the first file
                    Path(target_file).unlink(missing_ok=True)
                    Path(target_file).symlink_to(first_file)
                    # Write input file for raptor (space separated)
                    filehibf.write(target_file + " " + " ".join(files[1:]) + "\n")
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

            tx = time.time()
            print_log("raptor prepare", cfg.quiet)
            run_raptor_prepare_cmd = " ".join([cfg.path_exec['raptor'], "prepare",
                                               "--input '" + raptor_input_file + "'",
                                               "--output '" + build_output_folder + "'",
                                               "--kmer " + str(cfg.kmer_size),
                                               "--window " + str(cfg.window_size),
                                               "--quiet" if not cfg.verbose else "",
                                               "--threads " + str(cfg.threads),
                                               "--version-check 0"])
            run(run_raptor_prepare_cmd, quiet=cfg.quiet)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

            tx = time.time()
            print_log("raptor layout", cfg.quiet)
            # Use info file as input for raptor 
            run_raptor_layout_cmd = " ".join([cfg.path_exec['raptor'], "layout",
                                              "--input-file '" + build_output_folder + "minimiser.list'",
                                              "--tmax " + str(math.ceil(math.sqrt(n_bins) /64.0 ) * 64),
                                              "--num-hash-functions " + str(cfg.hash_functions),
                                              "--false-positive-rate " + str(cfg.max_fp),
                                              "--output-filename '" + build_output_folder + "raptor_layout.binning.out'",
                                              "--threads " + str(cfg.threads),
                                              "--version-check 0"])
            run(run_raptor_layout_cmd, quiet=cfg.quiet)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

            tx = time.time()
            print_log("raptor build", cfg.quiet)
            run_raptor_build_cmd = " ".join([cfg.path_exec['raptor'], "build",
                                             "--output '" + cfg.db_prefix + ".hibf" + "'",
                                             "--threads " + str(cfg.threads),
                                             "--quiet" if not cfg.verbose else "",
                                             "--input '" + build_output_folder + "raptor_layout.binning.out'",
                                             "--version-check 0"])
            run(run_raptor_build_cmd, quiet=cfg.quiet)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        else:
            tx = time.time()
            print_log("Building index (ganon-build)", cfg.quiet)
            run_ganon_build_cmd = " ".join([cfg.path_exec['build'],
                                            "--input-file '" + target_info_file + "'",
                                            "--output-file '" + cfg.db_prefix + ".ibf" + "'",
                                            "--kmer-size " + str(cfg.kmer_size),
                                            "--window-size " + str(cfg.window_size),
                                            "--hash-functions " + str(cfg.hash_functions),
                                            "--mode " + cfg.mode,
                                            "--min-length " + str(cfg.min_length),
                                            "--max-fp " + str(cfg.max_fp) if cfg.max_fp else "",
                                            "--filter-size " + str(cfg.filter_size) if cfg.filter_size else "",
                                            "--tmp-output-folder '" + build_output_folder + "'",
                                            "--threads " + str(cfg.threads),
                                            "--verbose" if cfg.verbose else "",
                                            "--quiet" if cfg.quiet else ""])
            run(run_ganon_build_cmd, quiet=cfg.quiet)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        save_state(which_call + "_run", files_output_folder)

    # Set output database files
    db_files_ext = []
    if cfg.filter_type == "hibf":
        db_files_ext.append("hibf")
    else:
        db_files_ext.append("ibf")
    if cfg.taxonomy != "skip":
        db_files_ext.append("tax")

    print_log("Database: " + ", ".join([cfg.db_prefix + "." + e for e in db_files_ext]), cfg.quiet)
    if all([check_file(cfg.db_prefix + "." + e) for e in db_files_ext]):
        # Do not delete temp (mostly tests, hidden param)
        if not cfg.keep_files:
            if which_call == "build_custom":
                # remove temporary files folder
                shutil.rmtree(files_output_folder)
            else:
                # remove tmp build folder
                shutil.rmtree(build_output_folder, ignore_errors=True)

        # remove save states
        clear_states(which_call, files_output_folder)
        print_log("Build finished successfully", cfg.quiet)
        return True
    else:
        print_log("ERROR: build failed - one or more database files not found or empty", cfg.quiet)
        return False


########################################################################################################################


def parse_input_file(input_file, info_cols, input_target, quiet):
    """
    parse user provided --input-file with all specifications for input
    """
    info = pd.read_csv(input_file,
                       sep="\t",
                       header=None,
                       skiprows=0,
                       dtype=object)
    # Need to rename columns after parsing since file can have 1 to 5
    info.rename(columns=lambda x: info_cols[x], inplace=True)
    info = pd.concat([info, pd.DataFrame(columns=info_cols)])
    
    # If no target was provided and target is file, use filename
    if info["target"].isna().all() and input_target == "file":
        info["target"] = info["file"].apply(os.path.basename)

    # Validate files provided and remove not found
    total_files = len(info["file"].unique().tolist())
    valid_files = validate_input_files(info["file"].unique().tolist(), "", quiet)
    if total_files - len(valid_files) > 0:
        info = info[info["file"].isin(valid_files)]
        print_log(" - " + str(total_files - len(valid_files)) + " invalid files skipped", quiet)
    return info


def load_input(cfg, input_files, build_output_folder):
    """
    Load basic target info, either provided as --target-info
    or extracted from file/sequences
    """
    tx = time.time()
    info_cols = ["file", "target", "node", "specialization", "specialization_name"]

    # Parse/load info without setting index
    if cfg.input_file:
        print_log("Parsing --input-file " + cfg.input_file, cfg.quiet)
        info = parse_input_file(cfg.input_file, info_cols, cfg.input_target, cfg.quiet)

        # If --input-target sequence, split files by sequence and validate headers
        if cfg.input_target == "sequence":
            print_log("\nMatching --input-file entries with sequences", cfg.quiet)
            shape_tmp = info.shape[0]
            # Split files and return headers obtained
            info_seqs = parse_sequence_accession(info["file"].unique().tolist(), info_cols, build_output_folder)
            
            # Merge with provided file, keeping only matching entries
            info = pd.merge(left=info, right=info_seqs, on="target", suffixes=("", "_seqs"))[info_cols + ["file_seqs"]]
            info["file"] = info["file_seqs"]
            info.drop("file_seqs", axis=1, inplace=True)

            if shape_tmp - info.shape[0] > 0:
                print_log(" - " + str(shape_tmp - info.shape[0]) + " invalid sequences skipped", cfg.quiet)

    else:
        
        if cfg.input_target == "sequence":
            print_log("Parsing sequences from --input (" + str(len(input_files)) + " files)", cfg.quiet)
            info = parse_sequence_accession(input_files, info_cols, build_output_folder)
        else:
            print_log("Parsing files from --input (" + str(len(input_files)) + " files)", cfg.quiet)
            info = parse_file_accession(input_files, info_cols)

    # Drop cols without values
    shape_tmp = info.shape[0]
    info.dropna(how="all", inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " invalid entries skipped", cfg.quiet)

    # Drop cols without target
    shape_tmp = info.shape[0]
    info.dropna(subset=["target"], inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " invalid targets skipped", cfg.quiet)

    # Drop duplicated target
    shape_tmp = info.shape[0]
    info.drop_duplicates(subset=['target'], inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " duplicated targets skipped", cfg.quiet)

    # set target as index
    info.set_index('target', inplace=True)
    print_log(" - " + str(info.shape[0]) + " unique entries", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return info


def load_taxonomy(cfg, build_output_folder):
    """
    load/download chosen taxonomy from multitax
    """
    tx = time.time()

    if cfg.taxonomy_files:
        print_log("Parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)
    else:
        print_log("Downloading and parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)

    if cfg.taxonomy == "ncbi":
        tax = NcbiTx(files=cfg.taxonomy_files)
    elif cfg.taxonomy == "gtdb":
        tax = GtdbTx(files=cfg.taxonomy_files, output_prefix=build_output_folder)

    # If level is not in special targets or leaves and present in available ranks
    if cfg.level not in [None, "leaves"] + cfg.choices_level:
        if cfg.level not in set(tax._ranks.values()):
            print_log(" - " + cfg.level + " not found in taxonomic ranks, changing to --level 'leaves'", cfg.quiet)
            cfg.level = 'leaves'

    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    return tax


def write_tax(tax_file, info, tax, genome_sizes, user_bins_col, level, input_target):
    """
    write tabular taxonomy file .tax
    may include specialization as nodes
    """

    # Add specialization if not using direct taxonomic nodes
    if user_bins_col != "node":
        # Set rank to level or input_target
        tax_rank = level if level else input_target
        for target, row in info.iterrows():
            tax_node = row["specialization"] if user_bins_col == "specialization" else target
            tax_name = row["specialization_name"] if user_bins_col == "specialization" else target

            # Check if node is already present with correct parent
            # in case of input-target sequence, info has repeated pairs of node/parent
            if tax.latest(tax_node) is tax.undefined_node:
                tax.add(tax_node, row["node"], name=tax_name, rank=tax_rank)
            else:
                assert tax.parent(tax_node)==row["node"]
            
    # Write filtered taxonomy with added nodes
    rm_files(tax_file)
    tax.write(tax_file)

    # add genome_sizes col, either from node or parent (specialization)
    tax_df = pd.read_csv(tax_file, names=["node", "parent", "rank", "name"], delimiter='\t', dtype=str)
    # Get estimated genome size from parent in case of specialization 
    tax_df["genome_size"] = tax_df.apply(lambda d: genome_sizes[d.node] if d.node in genome_sizes else genome_sizes[d.parent], axis=1)
    tax_df.to_csv(tax_file, sep="\t", header=False, index=False)

def write_target_info(info, user_bins_col, target_info_file):
    """
    write tabular file to be parsed by ganon-build with: file <tab> target
    """
    with open(target_info_file, "w") as outf:
        for target, row in info.iterrows():
            t = row[user_bins_col] if user_bins_col != "target" else target
            print(row["file"], t, sep="\t", end="\n", file=outf)


def write_info_file(info, filename):
    """
    write tabular file to be re-used as --input-file (sort cols in the right order)
    db_prefix.info.tsv
    """
    info.reset_index()[['file',
                        'target',
                        'node',
                        'specialization',
                        'specialization_name']].to_csv(filename,
                                                       sep="\t",
                                                       header=False,
                                                       index=False)


def validate_specialization(info, quiet):
    """
    validate specialization for each node
    each specialization can have only one parent node
    and invalid nodes
    """

    tx = time.time()
    print_log("Validating specialization", quiet)
    # if all entries are null, no specialization was retrieved
    if all(info.specialization.isna()):
        print_log(" - No specialization provided/retrieved", quiet)
    else:
        # check for invalid specialization entries
        idx_null_spec = info.specialization.isna()

        # get unique tuples node-specialization
        node_spec = info[['node', 'specialization']].drop_duplicates()

        # check for duplicated specialization in the tuples
        idx_multi_parent_spec = info.specialization.isin(node_spec.specialization[node_spec.specialization.duplicated(keep=False)].unique())

        # merge indices for invalid entries
        idx_replace = idx_null_spec | idx_multi_parent_spec

        if idx_replace.any():
            # replace invalid specialization entries with target
            info.loc[idx_replace, "specialization"] = info.index[idx_replace]
            info.loc[idx_replace, "specialization_name"] = info.index[idx_replace]
            print_log(" - " + str(sum(idx_replace)) + " invalid specialization entries replaced by target", quiet)

    # Skip invalid nodes
    shape_tmp = info.shape[0]
    info.dropna(subset=["specialization"], inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " entries without valid specialization skipped", quiet)

    # Fill names not provided with specialization
    info["specialization_name"].fillna(info["specialization"], inplace=True)

    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", quiet)


def validate_taxonomy(info, tax, cfg):
    """
    validate taxonomy: convert to latest nodes (tax.latest)
    and chosen level (tax.parent_rank)
    """
    tx = time.time()
    print_log("Validating taxonomy", cfg.quiet)

    # Get latest and valid taxonomic nodes
    info["node"] = info["node"].apply(tax.latest)

    # If level is set and not leaves or reserved
    if cfg.level and cfg.level not in ["leaves"] + cfg.choices_level:
        info["node"] = info["node"].apply(lambda n: tax.parent_rank(n, cfg.level))

    # Skip invalid nodes (na == tax.undefined_node (None))
    na_entries = info["node"].isna().sum()

    if na_entries > 0:
        info.dropna(subset=["node"], inplace=True)
        print_log(" - " + str(na_entries) + " entries without valid taxonomic nodes skipped", cfg.quiet)

    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)


def get_gu_current_version(assembly_summary):
    """
    return current version from genome_updater
    """
    return os.path.dirname(os.readlink(assembly_summary))


def restart_build(fld):
    """
    delete temporary folder to start build from scratch
    """
    if os.path.isdir(fld):
        shutil.rmtree(fld)


def restart_update(fld):
    """
    delete save states to start update from scratch
    """
    clear_states("update", fld)


def clear_states(prefix, folder):
    """
    delete all build/update save states
    """
    rm_files([folder + prefix + "_download",
              folder + prefix + "_parse",
              folder + prefix + "_run"])


def save_config(cfg, config_file):
    """
    save configuration for updates based on an instance of the Config class
    """
    v = vars(cfg)
    v["version"] = cfg.version
    with open(config_file, "wb") as file:
        pickle.dump(v, file)


def load_config(config_file):
    """
    load configuration
    """
    with open(config_file, "rb") as file:
        cfg = pickle.load(file) 
    return cfg
