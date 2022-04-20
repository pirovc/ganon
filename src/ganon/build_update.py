import time
import pandas as pd
import os

from ganon.config import Config
from ganon.util import print_log
from ganon.util import set_out_folder
from ganon.util import set_out_files
from ganon.util import run
from ganon.util import validate_input_files
from ganon.util import rm_folder
from ganon.tax_util import get_file_info
from ganon.tax_util import get_sequence_info
from ganon.tax_util import parse_sequence_accession
from ganon.tax_util import parse_file_accession

from multitax import NcbiTx, GtdbTx


def build(cfg):
    gu_output_folder = cfg.db_prefix + "_files/"
    resume_download = False

    if os.path.exists(gu_output_folder) and not cfg.restart:
        print_log("Download folder found [" + gu_output_folder + "], trying to resume\n", cfg.quiet)
        resume_download = True
    elif not set_out_folder(gu_output_folder, cfg.restart):
        return False

    tx = time.time()
    print_log("Downloading files from " + cfg.source + " [" + ",".join(cfg.organism_group) + "]", cfg.quiet)
    run_genome_updater_cmd = " ".join([cfg.path_exec['genome_updater'],
                                       "-d " + cfg.source,
                                       "-g " + ",".join(cfg.organism_group),
                                       "-A " + str(cfg.top) if cfg.taxonomy == "ncbi" else "",
                                       "-l 'complete genome'" if cfg.complete_genomes else "",
                                       "-f genomic.fna.gz",
                                       "-t " + str(cfg.threads),
                                       "-o " + gu_output_folder,
                                       "-z " if cfg.taxonomy == "gtdb" else "",
                                       "-m",
                                       "-i" if resume_download else "",
                                       "-s" if cfg.quiet else "-w",
                                       cfg.genome_updater if cfg.genome_updater else ""])
    run(run_genome_updater_cmd, quiet=cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # get current version from assembly_summary
    assembly_summary = gu_output_folder + "assembly_summary.txt"
    input_folder = gu_output_folder + os.path.dirname(os.readlink(assembly_summary)) + "/files/"

    build_custom_params = {"input": input_folder,
                           "input_extension": "fna.gz",
                           "input_target": "file",
                           "level": "name",
                           "ncbi_file_info": assembly_summary}

    build_default_params = {"db_prefix": cfg.db_prefix,
                            "taxonomy": cfg.taxonomy,
                            "threads": cfg.threads,
                            "max_fp": cfg.max_fp,
                            "filter_size": cfg.filter_size,
                            "kmer_size": cfg.kmer_size,
                            "window_size": cfg.window_size,
                            "hash_functions": cfg.hash_functions,
                            "restart": cfg.restart,
                            "verbose": cfg.verbose,
                            "quiet": cfg.quiet,
                            "ganon_path": cfg.ganon_path,
                            "n_refs": cfg.n_refs,
                            "n_batches": cfg.n_batches}

    build_custom_params.update(build_default_params)
    return build_custom(cfg=Config("build-custom", **build_custom_params))


def update(cfg):
    return True


def build_custom(cfg):

    tax = None
    input_files = []
    tmp_output_folder = cfg.db_prefix + "_tmp/"

    # Define if taxonomy/specialization is used as a target to build the database
    use_spec_target = True if cfg.level in cfg.choices_level else False

    # Set working folder, check out files
    if not set_out_files(cfg.db_prefix, ["ibf", "tax"], cfg.restart):
        return False
    if not set_out_folder(tmp_output_folder, cfg.restart):
        return False

    # Retrieve and check input files or folders
    input_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet)
    if not input_files:
        print_log("ERROR: No valid input files found")
        return False
    print_log("")

    # Set --input-target if not manually set
    if not cfg.input_target:
        cfg.input_target = "sequence" if len(input_files) == 1 else "file"

    # Set-up taxonomy
    if cfg.taxonomy != "none":
        tx = time.time()

        if cfg.taxonomy_files:
            print_log("Parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)
        else:
            print_log("Downloading and parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)

        if cfg.taxonomy == "ncbi":
            tax = NcbiTx(files=cfg.taxonomy_files)
        elif cfg.taxonomy == "gtdb":
            tax = GtdbTx(files=cfg.taxonomy_files, output_prefix=tmp_output_folder + "gtdb_")

        # If level is not in special targets or leaves and present in available ranks
        if cfg.level not in [None, "leaves"] + cfg.choices_level:
            if cfg.level not in set(tax._ranks.values()):
                print_log(" - " + cfg.level + " not found in taxonomic ranks, changing to --level 'leaves'", cfg.quiet)
                cfg.level = 'leaves'

        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # Set-up input info
    info = load_input(cfg, input_files)

    # Retrieve target info if taxonomy or specialization is required (and if file is not provided)
    if (tax or use_spec_target) and not cfg.input_file:
        if cfg.input_target == "sequence":
            get_sequence_info(cfg, info, tax, tmp_output_folder, use_spec_target)
        else:
            get_file_info(cfg, info, tax, tmp_output_folder)

    # Validate taxonomic node only if taxonomy is provided
    if tax:
        validate_taxonomy(info, tax, cfg)
        if info.empty:
            print_log("ERROR: Could not retrieve any taxonomic information")
            return False

    # Validate specialization (required after taxonomy)
    if use_spec_target:
        validate_specialization(info, cfg.quiet)
        if info.empty:
            print_log("ERROR: Could not retrieve any specialization information")
            return False

    print(info)

    # Define target fro user bins
    user_bins = "target"
    if use_spec_target:
        user_bins = "specialization"
    elif cfg.level:  # if any level is provided
        user_bins = "node"

    # Filter and write taxonomy
    if tax:
        # filter only used tax. nodes
        tax.filter(info["node"].unique())
        write_tax(cfg, info, tax, user_bins)

    # Write aux file for ganon-build
    write_target_info(info, cfg.input_target, user_bins, tmp_output_folder)

    # run ganon-build

    #rm_folder(tmp_output_folder)
    return True


def load_input(cfg, input_files):
    """
    Load basic target info, either provided as --target-info
    or extracted from file/sequences
    """
    target_info_columns = ["target", "node", "specialization", "file"]
    tx = time.time()
    if cfg.input_file:
        print_log("Parsing --target-info " + cfg.input_file, cfg.quiet)
        info = pd.read_csv(cfg.input_file,
                           sep="\t",
                           header=None,
                           skiprows=0,
                           index_col=None,
                           dtype="str",
                           names=target_info_columns)
    else:
        if cfg.input_target == "sequence":
            print_log("Parsing sequences from --input (" + str(len(input_files)) + " files)", cfg.quiet)
            info = parse_sequence_accession(input_files, target_info_columns)
        else:
            print_log("Parsing --input (" + str(len(input_files)) + " files)", cfg.quiet)
            info = parse_file_accession(input_files, target_info_columns)

    # Drop cols without values
    shape_tmp = info.shape[0]
    info.dropna(how="all", inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " invalid entries skipped", cfg.quiet)

    # Drop duplicated target (seqid or fileid)
    shape_tmp = info.shape[0]
    info.drop_duplicates(subset=['target'], inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " duplicated entries skipped", cfg.quiet)

    # set target as index
    info.set_index('target', inplace=True)
    print_log(" - " + str(info.shape[0]) + " unique entries", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return info


def write_tax(cfg, info, tax, user_bins):
    """
    write tabular taxonomy file .tax
    may include specialization as nodes
    """

    tax_file = cfg.db_prefix + ".tax"

    # Write filtered "standard" taxonomy
    tax.write(tax_file)

    # Set rank to level or input_target
    rank = cfg.level if cfg.level else cfg.input_target

    # Add specialization if not using direct taxonomic nodes
    if user_bins != "node":
        with open(tax_file, "a") as outf:
            for target, row in info.iterrows():
                t = row[user_bins] if user_bins != "target" else target
                print(t, row["node"], rank, t, sep="\t", end="\n", file=outf)


def write_target_info(info, input_target, user_bins, tmp_output_folder):
    """
    write tabular file to be parsed by ganon-build with info about file, target [and sequence]
    """

    target_info_file = tmp_output_folder + "target_info.tsv"

    # Add specilization nodes
    with open(target_info_file, "w") as outf:
        for target, row in info.iterrows():
            t = row[user_bins] if user_bins != "target" else target
            s = target if input_target == "sequence" else ""
            print(row["file"], t, s, sep="\t", end="\n", file=outf)


def validate_specialization(info, quiet):
    """
    validate specialization in terms of relation to
    taxonomy (each specialization can have only one parent node)
    and invalid nodes
    """

    tx = time.time()
    print_log("Validating level specialization", quiet)
    # if all entries are null, no specialization was retrieved
    if all(info.specialization.isnull()):
        print_log(" - No valid level specialization could be defined\n", quiet)
    else:
        # check for invalid specialization entries
        idx_null_spec = info.specialization.isnull()

        # get unique tuples node-specialization
        node_spec = info[['node', 'specialization']].drop_duplicates()

        # check for duplicated specialization in the tuples
        idx_multi_parent_spec = info.specialization.isin(node_spec.specialization[node_spec.specialization.duplicated(keep=False)].unique())

        # merge indices for invalid entries
        idx_replace = idx_null_spec | idx_multi_parent_spec

        if idx_replace.any():
            # replace invalid specialization entries with target
            info.loc[idx_replace, "specialization"] = info.index[idx_replace]
            print_log(str(sum(idx_replace)) + " invalid specialization replaced by target\n", quiet)

        # Skip invalid nodes (na == tax.undefined_node (None))
        shape_tmp = info.shape[0]
        info.dropna(subset=["specialization"], inplace=True)
        if shape_tmp - info.shape[0] > 0:
            print_log(" - " + str(shape_tmp - info.shape[0]) + " entries without valid specialization skipped", quiet)

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
    shape_tmp = info.shape[0]
    info.dropna(subset=["node"], inplace=True)
    if shape_tmp - info.shape[0] > 0:
        print_log(" - " + str(shape_tmp - info.shape[0]) + " entries without valid taxonomic nodes skipped", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)


def update_custom(cfg):
    return True
