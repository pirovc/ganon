import time, math
import pandas as pd
import numpy as np
import taxsbp.taxsbp
from io import StringIO
from ganon.bins import Bins
from ganon.gnn import Gnn
from ganon.seqinfo import SeqInfo
from ganon.tax import Tax
from ganon.util import *


def build(cfg):
    # validate input files
    input_files = validate_input_files(cfg.input_files, cfg.input_directory, cfg.input_extension, cfg.quiet)
    if len(input_files)==0:
        print_log("ERROR: No valid input files found", cfg.quiet)
        return False

    # Set db prefixes
    db_prefix = {prefix: cfg.db_prefix + "." + prefix for prefix in  ["ibf", "map", "tax", "gnn"]}

    # Set temporary working folder 
    tmp_output_folder = cfg.db_prefix + "_tmp/"
    if not set_tmp_folder(tmp_output_folder): return False

    # Set up taxonomy
    ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(cfg.taxdump_file, tmp_output_folder, cfg.quiet)

    # Parse .tax    
    tx = time.time()
    print_log("Parsing taxonomy", cfg.quiet)
    tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # load seqinfo (file or seqids)
    seqinfo = load_seqinfo(cfg, input_files)

    # Retrieve sequence information
    if not cfg.seq_info_file:
        retrieve_seqinfo(seqinfo, tmp_output_folder, input_files, cfg)

    # Check for valid specialization
    if cfg.specialization:
        replaced_spec = seqinfo.validate_specialization()
        if replaced_spec:
            print_log(str(replaced_spec) + " invalid specialization entries (sequence accession used instead)\n", cfg.quiet)

    # Write seq-info-file
    if not cfg.seq_info_file and cfg.write_seq_info_file:
        seqinfo.write(cfg.db_prefix + ".seqinfo.txt")

    # check sequences compared to bins
    added_seqids, _, _ = check_updated_seqids(set(seqinfo.get_seqids()), set())
    # Ignore removed sequences if not doing complete update
    print_log("Build: adding " + str(len(added_seqids)) + " sequences", cfg.quiet)
    print_log("", cfg.quiet)

    if not added_seqids:
        print_log("No valid seq. info to build", cfg.quiet)
        rm_tmp_folder(tmp_output_folder)
        return False

    # Set or calculate best --bin-length
    if cfg.bin_length:
        bin_length = cfg.bin_length
    else:
        tx = time.time()
        print_log("Simulating parameters", cfg.quiet)
        bin_length = estimate_bin_length(cfg, seqinfo, tax)
        if cfg.verbose:
            print_log(" - --bin-length " + str(bin_length), cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # Set fragment length
    if cfg.fragment_length == -1:  # if ==-1 set default
        fragment_length = bin_length - cfg.overlap_length
    elif cfg.fragment_length == 0:  # if ==0 deactivate
        fragment_length = 0
    else:   # user input
        fragment_length = cfg.fragment_length - cfg.overlap_length

    tx = time.time()
    print_log("Running taxonomic clustering (TaxSBP)", cfg.quiet)
    bins = run_taxsbp(seqinfo, bin_length, fragment_length, cfg.overlap_length, cfg.rank, cfg.specialization, ncbi_nodes_file, ncbi_merged_file, cfg.verbose)
    # bin statistics
    actual_number_of_bins = bins.get_number_of_bins()
    optimal_number_of_bins = optimal_bins(actual_number_of_bins)
    max_length_bin = bins.get_max_bin_length()
    max_kmer_count = estimate_elements(max_length_bin, cfg.kmer_size, cfg.window_size)
    max_split_bins = bins.get_max_split_bins(use_specialization=True if cfg.specialization else False)
    if cfg.verbose:
        print_log(" - " + str(actual_number_of_bins) + " bins created", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # Get approximate optimal parameters from user input and taxsbp result
    optimal_params = derive_bf_params(max_kmer_count, cfg.max_fp, 0, cfg.hash_functions)

    # When fixed size is too small
    if optimal_params["hash_functions"] == 0:
        optimal_params["hash_functions"] = 3

    if cfg.verbose:
        print_log("IBF params:", cfg.quiet)
        try:
            correction_ratio = split_correction_ratio(optimal_params["false_positive"], optimal_params["hash_functions"], max_split_bins)
        except:
            correction_ratio = 1

        print_log(" - # hash functions: " + str(optimal_params["hash_functions"]), cfg.quiet)
        print_log(" - # targets: " + str(len(bins.get_specialization() if cfg.specialization else bins.get_taxids())), cfg.quiet)
        print_log(" - # bins/optimal bins: " + str(actual_number_of_bins) + "/" + str(optimal_number_of_bins), cfg.quiet)
        print_log(" - target with most splits/bin: " + str(max_split_bins), cfg.quiet)
        print_log(" - max. correction ratio: " + str("{0:.5f}".format(correction_ratio)) , cfg.quiet)
        if cfg.filter_size:
            print_log(" - fixed filter size: " + str("{0:.2f}".format(cfg.filter_size)) + "MB" , cfg.quiet)
        else:
            print_log(" - max. false positive: " + str("{0:.5f}".format(cfg.max_fp)), cfg.quiet)
            print_log(" - upper limit filter size: " + str("{0:.2f}".format(bits2mb(optimal_params["size_bits"] * optimal_number_of_bins * correction_ratio))) + "MB", cfg.quiet)
            if(cfg.window_size>0):
                print_log(" ** real filter size should be considerably smaller", cfg.quiet)
        print_log("")

    # Build database files (map, tax, gnn)
    tx = time.time()
    print_log("Building database files", cfg.quiet)
    # Write .map file
    print_log(" - " + db_prefix["map"], cfg.quiet)
    bins.write_map_file(db_prefix["map"], use_specialization=True if cfg.specialization else False)

    # Write .tax file
    print_log(" - " + db_prefix["tax"], cfg.quiet)
    # filter only used taxids
    tax.filter(bins.get_taxids())
    # add specialization nodes
    if cfg.specialization:
        tax.add_nodes(bins.get_specialization_taxid(), cfg.specialization)
    tax.write(db_prefix["tax"])

    if cfg.specialization and cfg.rank != "leaves":
        cfg.rank = "leaves"
        if cfg.verbose:
            print_log(" - --rank is set to leaves when using specialization values", cfg.quiet)
        
    # Write .gnn file
    print_log(" - " + db_prefix["gnn"], cfg.quiet)
    gnn = Gnn(kmer_size=cfg.kmer_size,
              window_size=cfg.window_size,
              hash_functions=optimal_params["hash_functions"],
              number_of_bins=actual_number_of_bins,
              rank=cfg.rank,
              specialization=cfg.specialization,
              bin_length=bin_length,
              fragment_length=fragment_length,
              overlap_length=cfg.overlap_length,
              bins=bins.get_list())
    gnn.write(db_prefix["gnn"])
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    print_log("Building index (ganon-build)", cfg.quiet)
    # Write aux. file for ganon
    acc_bin_file = tmp_output_folder + "acc_bin.txt"
    bins.write_acc_bin_file(acc_bin_file)

    # Free memory for build
    del seqinfo
    del bins
    del tax
    del gnn

    run_ganon_build_cmd = " ".join([cfg.path_exec['build'],
                                    "--seqid-bin-file " + acc_bin_file,
                                    "--filter-size-mb " + str(cfg.filter_size) if cfg.filter_size else "--false-positive " + str(cfg.max_fp),
                                    "--kmer-size " + str(cfg.kmer_size),
                                    "--window-size " + str(cfg.window_size) if cfg.window_size else "",
                                    "--count-hashes " if cfg.window_size else "",
                                    "--map " + db_prefix["map"],
                                    "--hash-functions " + str(optimal_params["hash_functions"]),
                                    "--threads " + str(cfg.threads),
                                    "--output-filter-file " + db_prefix["ibf"],
                                    "--verbose" if cfg.verbose else "",
                                    "--quiet" if cfg.quiet else "",
                                    "--n-refs " + str(cfg.n_refs) if cfg.n_refs is not None else "",
                                    "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                    "--reference-files " + ",".join([file for file in input_files]) if input_files and not cfg.input_directory else "",
                                    "--directory-reference-files " + cfg.input_directory if cfg.input_directory else "",
                                    "--extension " + cfg.input_extension if cfg.input_extension else ""])
    stdout, stderr = run(run_ganon_build_cmd, print_stderr=True)

    # Delete temp files
    rm_tmp_folder(tmp_output_folder)

    return True


def update(cfg):
    tx = time.time()

    # validate input files
    input_files = validate_input_files(cfg.input_files, cfg.input_directory, cfg.input_extension, cfg.quiet)
    if len(input_files)==0:
        print_log("ERROR: No valid input files found", cfg.quiet)
        return False

    # Set db prefixes
    db_prefix = {prefix:cfg.db_prefix + "." + prefix for prefix in  ["ibf","map","tax","gnn"]}  

    # Set temporary working folder (current or new output)
    tmp_output_folder = cfg.output_db_prefix + "_tmp/" if cfg.output_db_prefix else cfg.db_prefix + "_tmp/"
    if not set_tmp_folder(tmp_output_folder): return False

    # Load .gnn file   
    gnn = Gnn(file=db_prefix["gnn"])

    # If specialization was set on database
    if gnn.specialization:
        # if not provided by user, use defition of database
        if not cfg.specialization: cfg.specialization=gnn.specialization
        print_log("Using --specialization " + cfg.specialization, cfg.quiet)
    else:
        if cfg.specialization: 
            # If user defined specialization on update but database has none
            print_log("ERROR: not possible to update a database with --specialization if it was built without it", cfg.quiet)
            return False

    # load bins
    bins = Bins(taxsbp_ret=gnn.bins, use_specialization=True if cfg.specialization else False)

    # load seqinfo (file or seqids)
    seqinfo = load_seqinfo(cfg, input_files)

    # check sequences compared to bins
    added_seqids, removed_seqids, kept_seqids = check_updated_seqids(set(seqinfo.get_seqids()), set(bins.get_seqids()))
    # Ignore removed sequences if not doing complete update
    if cfg.update_complete: 
        print_log("Update: adding " + str(len(added_seqids)) + " sequences, removing " + str(len(removed_seqids)) + " sequences, keeping " + str(len(kept_seqids)) + " sequences", cfg.quiet)
    else:
        removed_seqids = []
        print_log("Update: adding " + str(len(added_seqids)) + " sequences, ignoring " + str(len(kept_seqids)) + " repeated sequences", cfg.quiet)
    print_log("", cfg.quiet)

    if not added_seqids and not removed_seqids:
        print_log("Nothing to update", cfg.quiet)
        rm_tmp_folder(tmp_output_folder)
        return False

    if cfg.update_complete:
        # Remove already included seqids to just retrieve information for added sequences
        seqinfo.remove_seqids(kept_seqids | removed_seqids)
    else:
        # Remove seqids already present in the current version (repeated entries)
        seqinfo.remove_seqids(kept_seqids)

    # retrive sequence information (after removing invalid seqids)
    if not cfg.seq_info_file:
        retrieve_seqinfo(seqinfo, tmp_output_folder, input_files, cfg)

    # Convert cols data types
    if cfg.specialization:
        replaced_spec = seqinfo.validate_specialization()
        if replaced_spec:
            print_log(str(replaced_spec) + " invalid specialization entries (sequence accession used instead)\n", cfg.quiet)

    if not cfg.seq_info_file and cfg.write_seq_info_file:
        seqinfo.write(cfg.output_db_prefix+".seqinfo.txt")

    # save set of current binids
    previous_binids = set(bins.get_binids())
    # remove seqids from bins if performing update complete
    if cfg.update_complete and removed_seqids:
        bins.remove_seqids(removed_seqids)
    # save set of kept binids after removal
    kept_binids = set(bins.get_binids())

    # Set up taxonomy files
    ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(cfg.taxdump_file, tmp_output_folder, cfg.quiet)

    tx = time.time()
    print_log("Running taxonomic clustering (TaxSBP)", cfg.quiet)
    updated_bins = run_taxsbp(seqinfo, gnn.bin_length, gnn.fragment_length, gnn.overlap_length, gnn.rank, cfg.specialization, ncbi_nodes_file, ncbi_merged_file, cfg.verbose, bins=bins)
    # bin statistics
    taxsbp_binids = set(updated_bins.get_binids())
    removed_binids = previous_binids.difference(kept_binids | taxsbp_binids)
    new_binids = taxsbp_binids.difference(previous_binids)
    updated_binids = taxsbp_binids.intersection(previous_binids)
    print_log(" - " + str(len(new_binids)) + " bins added, " + str(len(updated_binids)) + " bins updated, " + str(len(removed_binids)) + " bins removed", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    tx = time.time()
    print_log("Updating database files", cfg.quiet)
    # load new taxonomy
    print_log(" - " + cfg.output_db_prefix + ".tax" if cfg.output_db_prefix else db_prefix["tax"], cfg.quiet)
    tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
    # Update and write .tax file
    # filter only used taxids
    tax.filter(updated_bins.get_taxids())
    # add specialization nodes
    if cfg.specialization:
        tax.add_nodes(updated_bins.get_specialization_taxid(), cfg.specialization)

    # Load old .tax file into new taxonomy
    tax.merge(Tax([db_prefix["tax"]]))
    # Write .tax file
    tax.write(cfg.output_db_prefix + ".tax" if cfg.output_db_prefix else db_prefix["tax"])
    # TODO - remove entries from .tax from removed entries of the db

    # merge updated and old bins together
    bins.merge(updated_bins)

    # Write .gnn file
    print_log(" - " + cfg.output_db_prefix + ".gnn" if cfg.output_db_prefix else db_prefix["gnn"], cfg.quiet)
    gnn.bins = bins.get_list()  # save updated bins
    gnn.number_of_bins = bins.get_number_of_bins()  # add new bins count
    # set new specialization to gnn
    gnn.specialization = cfg.specialization
    gnn.write(cfg.output_db_prefix + ".gnn" if cfg.output_db_prefix else db_prefix["gnn"])

    # Recreate .map file based on the new bins
    print_log(" - " + cfg.output_db_prefix + ".map" if cfg.output_db_prefix else db_prefix["map"], cfg.quiet)
    bins.write_map_file(cfg.output_db_prefix + ".map" if cfg.output_db_prefix else db_prefix["map"], use_specialization=True if cfg.specialization else False)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    tx = time.time()
    print_log("Updating index (ganon-build)", cfg.quiet)

    # Write aux. file for ganon
    # This file has to contain all new sequences
    # in case of update_complete
    acc_bin_file = tmp_output_folder + "acc_bin.txt"

    if cfg.update_complete:
        # all sequences from the bins with added/removed sequences should be written
        bins.write_acc_bin_file(acc_bin_file, new_binids | updated_binids)
        # If all sequences of a bin were removed and no new sequence added
        # insert a dummy entry for ganon-build to clear the bin
        if removed_binids:
            with open(acc_bin_file, "a") as abf:
                for b in removed_binids:
                    print(0, 0, 0, b, sep="\t", file=abf)
    else:
        # Only new sequences (updated_bins) either on old or new binids
        updated_bins.write_acc_bin_file(acc_bin_file)

    # Update with same values used for build
    kmer_size = gnn.kmer_size
    window_size = gnn.window_size
    hash_functions = gnn.hash_functions

    # Free memory for build
    del seqinfo
    del bins
    del updated_bins
    del tax
    del gnn

    # Temporary output filter 
    tmp_db_prefix_ibf = tmp_output_folder + "ganon.ibf"
    run_ganon_build_cmd = " ".join([cfg.path_exec['build'],
                                    "--update-filter-file " + db_prefix["ibf"],
                                    "--kmer-size " + str(kmer_size),
                                    "--window-size " + str(window_size) if window_size else "",
                                    "--count-hashes " if window_size else "",
                                    "--hash-functions " + str(hash_functions),
                                    "--seqid-bin-file " + acc_bin_file,
                                    "--output-filter-file " + tmp_db_prefix_ibf,
                                    "--threads " + str(cfg.threads),
                                    "--verbose" if cfg.verbose else "",
                                    "--quiet" if cfg.quiet else "",
                                    "--n-refs " + str(cfg.n_refs) if cfg.n_refs is not None else "",
                                    "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                    "--reference-files " + ",".join([file for file in input_files]) if input_files and not cfg.input_directory else "",
                                    "--directory-reference-files " + cfg.input_directory if cfg.input_directory else "",
                                    "--extension " + cfg.input_extension if cfg.input_extension else "",
                                    "--update-complete" if cfg.update_complete else ""])
    stdout, stderr = run(run_ganon_build_cmd, print_stderr=True)

    # move IBF to final location
    shutil.move(tmp_db_prefix_ibf, cfg.output_db_prefix + ".ibf" if cfg.output_db_prefix else db_prefix["ibf"])

    # Delete temp files
    rm_tmp_folder(tmp_output_folder)

    return True

def check_updated_seqids(new_seqids, old_seqids):
    # remove repeated from old bins
    added_seqids = new_seqids.difference(old_seqids)
    removed_seqids = old_seqids.difference(new_seqids)
    kept_seqids = old_seqids.difference(removed_seqids)

    return added_seqids, removed_seqids, kept_seqids

def parse_seqids(seqinfo, input_files, specialization, quiet, get_length: bool):
    for file in input_files:
        if get_length:
            # cat | zcat | gawk -> compability with osx
            run_get = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{if (seqlen){{print seqlen}}; printf substr($1,2)\"\\t\";seqlen=0;next;}} {{seqlen+=length($0)}}END{{print seqlen}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get, shell=True)
            parsed_stdout = pd.read_csv(StringIO(stdout), sep="\t", header=None, names=['seqid', 'length'])
        else:
            # cat | zcat | gawk -> compability with osx
            run_get = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{print substr($1,2)}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get, shell=True)
            parsed_stdout = pd.read_csv(StringIO(stdout), header=None, names=['seqid'])
        if specialization=="file":
            parsed_stdout["specialization"] = os.path.basename(file)
        elif specialization=="sequence":
            parsed_stdout["specialization"] = parsed_stdout["seqid"]
        seqinfo.append(parsed_stdout)

    # Drop duplicated seqids
    parsed_size = seqinfo.size()
    seqinfo.drop_duplicates()
    if seqinfo.size() < parsed_size: 
        print_log(" - " + str(parsed_size-seqinfo.size()) + " duplicated accessions were skipped", quiet)

    # Drop rows with zero length
    if get_length:
        parsed_size = seqinfo.size()
        seqinfo.drop_zeros(col='length')
        if seqinfo.size() < parsed_size: 
            print_log(" - " + str(parsed_size-seqinfo.size()) + " entries without sequence length were skipped", quiet)

def parse_eutils(seqinfo, tmp_output_folder, path_exec_get_seq_info, quiet, skip_len_taxid=False, get_assembly=False):
    seqid_file = tmp_output_folder + "seqids.txt"
    seqinfo.write_seqid_file(seqid_file)

    # always return all entries in the same order (-k)
    # report sequence accession if specialization not found (-r)
    run_get_seq_info_cmd = '{0} -k -r -i {1} {2} {3}'.format(path_exec_get_seq_info, seqid_file, "-a" if get_assembly else "", "-s" if skip_len_taxid else "")
    stdout, stderr = run(run_get_seq_info_cmd, print_stderr=True if not quiet else False, exit_on_error=False)

    # set "na" as NaN with na_values="na"
    if get_assembly:
        if skip_len_taxid:
            # Return only assembly
            seqinfo.paste_cols("specialization", pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=['seqid','assembly'], na_values="na")['assembly'])
        else:
            # Return full seqinfo
            seqinfo.clear()
            seqinfo.append(pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=['seqid','length','taxid','specialization'], na_values="na", dtype={'taxid': 'str'}))
    else:
        # Return seqid, len and taxid - keep specialization in the seqinfo
        seqinfo.paste_cols(["seqid","length","taxid"], pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=['seqid','length','taxid'], na_values="na", dtype={'taxid': 'str'}))

    # drop failed entries for lenght or taxid
    parsed_size = seqinfo.size()
    seqinfo.dropna(subset=["length","taxid"])
    if seqinfo.size() < parsed_size: 
        print_log(" - " + str(parsed_size-seqinfo.size()) + " entries could not be retrieved from eutils and were skipped", quiet)


def load_seqinfo(cfg, input_files):
    seqinfo = SeqInfo()
    if cfg.seq_info_file:
        tx = time.time()
        print_log("Parsing --seq-info-file", cfg.quiet)
        seqinfo.parse_seq_info_file(cfg.seq_info_file, use_specialization=True if cfg.specialization=="custom" else False)
        print_log(" - " + str(seqinfo.size()) + " unique sequence entries in the --seq-info-file " + cfg.seq_info_file, cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        tx = time.time()
        print_log("Extracting sequence identifiers", cfg.quiet)
        parse_seqids(seqinfo, input_files, cfg.specialization, cfg.quiet, get_length=False)
        print_log(" - " + str(seqinfo.size()) + " unique sequence headers successfully retrieved from " + str(len(input_files)) + " input file(s)" , cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return seqinfo


def retrieve_seqinfo(seqinfo, tmp_output_folder, input_files, cfg):
    # Max. # of sequences to use eutils as auto mode
    max_seqs_eutils = 50000
    # initialize total count
    seqid_total_count = seqinfo.size()

    # Define method to use
    if "auto" in cfg.seq_info_mode:
        if seqid_total_count > max_seqs_eutils: 
            seq_info_mode = ["nucl_gb", "nucl_wgs"]
        else:
            seq_info_mode = ["eutils"]
    elif "eutils" in cfg.seq_info_mode:
        seq_info_mode = ["eutils"]
    else:
        seq_info_mode = cfg.seq_info_mode

    if seq_info_mode[0]=="eutils":
        tx = time.time()
        print_log("Retrieving sequence information from NCBI E-utils", cfg.quiet)
        parse_eutils(seqinfo, tmp_output_folder, cfg.path_exec['get_seq_info'], cfg.quiet, skip_len_taxid=False, get_assembly=True if cfg.specialization=="assembly" else False)
        print_log(" - " + str(seqinfo.size()) + " sequences successfully retrieved", cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        # Retrieve seq. ids and lengths
        tx = time.time()
        print_log("Extracting sequence lengths", cfg.quiet)
        seqinfo.clear() # Clear seqinfo (get seqids again with length)
        parse_seqids(seqinfo, input_files, cfg.specialization, cfg.quiet, get_length=True)
        print_log(" - " + str(seqinfo.size()) + " sequences lenghts successfully retrieved", cfg.quiet)
        # Check if retrieved lengths are the same as number of inputs, reset counter
        if seqinfo.size() < seqid_total_count:
            print_log(" - could not retrieve lenght for " + str(seqid_total_count - seqinfo.size()) + " sequences", cfg.quiet)
            seqid_total_count = seqinfo.size()
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        dowloaded_acc2txid_files = []
        for acc2txid in seq_info_mode:
            dowloaded_acc2txid_files.append(get_accession2taxid(acc2txid, tmp_output_folder, cfg.quiet))

        print_log("Parsing accession2taxid files", cfg.quiet)
        count_acc2txid = parse_acc2txid(seqinfo, dowloaded_acc2txid_files)
        for acc2txid_file, cnt in count_acc2txid.items():
            print_log(" - " + str(cnt) + " entries found in the " + acc2txid_file.split("/")[-1] + " file", cfg.quiet)

        # filter out taxids not found
        seqinfo.dropna(subset=['taxid'])

        # Check if retrieved taxids are the same as number of inputs, reset counter
        if seqinfo.size() < seqid_total_count:
            print_log(" - could not retrieve taxid for " + str(seqid_total_count - seqinfo.size()) + " accessions", cfg.quiet)
            seqid_total_count = seqinfo.size()
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        if cfg.specialization=="assembly":
            tx = time.time()
            print_log("Retrieving assembly information from NCBI E-utils", cfg.quiet)
            parse_eutils(seqinfo, tmp_output_folder, cfg.path_exec['get_seq_info'], cfg.quiet, skip_len_taxid=True, get_assembly=True)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)


def get_accession2taxid(acc2txid, tmp_output_folder, quiet):
    tx = time.time()
    acc2txid_file = acc2txid + ".accession2taxid.gz"
    print_log("Downloading " + acc2txid_file, quiet)
    acc2txid_file = tmp_output_folder + acc2txid_file
    run_wget_acc2txid_file_cmd = 'wget -qO {0} "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/{1}.accession2taxid.gz"'.format(acc2txid_file, acc2txid)
    stdout, stderr = run(run_wget_acc2txid_file_cmd, print_stderr=True)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", quiet)
    return acc2txid_file


def parse_acc2txid(seqinfo, acc2txid_files):
    count_acc2txid = {}
    unique_seqids = set(seqinfo.get_seqids())
    for acc2txid in acc2txid_files:
        tmp_seqid_taxids = pd.read_csv(acc2txid,
                                       sep='\t',
                                       header=None,
                                       skiprows=1,
                                       usecols=[1, 2],
                                       names=['seqid', 'taxid'],
                                       index_col='seqid',
                                       converters={'seqid': lambda x: x if x in unique_seqids else None},
                                       dtype={'taxid': 'str'})
        tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids.index.notnull()] # keep only seqids used
        tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['taxid'] != "0"]  # filter out taxid==0

        # save count to return
        count_acc2txid[acc2txid] = tmp_seqid_taxids.shape[0]
        # merge taxid retrieved based on seqid (if any)
        if count_acc2txid[acc2txid]:
            seqinfo.join(tmp_seqid_taxids, "taxid")
        del tmp_seqid_taxids
        #if already found all seqids no need to parse all files till the end)
        if sum(count_acc2txid.values()) == len(unique_seqids):
            break

    return count_acc2txid


def bits2mb(v):
    return v / 8388608


def mb2bits(v):
    return v * 8388608


def optimal_bins(n):
    #return optimal number of bins for the IBF (multiples of 64)
    return math.ceil(n / float(64)) * 64


def estimate_n_bins(bin_len, overlap_len, groups_len):
    # Estimate an approximate number of bins give the parameters
    if bin_len > overlap_len:
        bins=[]
        for l in groups_len.values():
            bins.append(math.ceil(l/(bin_len-overlap_len)) if l > bin_len else 1)
        return sum(bins), max(bins)
    else:
        return 0, 0


def estimate_elements(bin_len, kmer_size, window_size):
    """
    Return estimation of elements given a bin lenght
    """
    if window_size > 0:
        # Estimate elements to be the median of min. and max. possible minimizers
        return bin_len - window_size + 1
    else:
        # Estimate elements to be all k-mers
        return bin_len - kmer_size + 1


def split_correction_ratio(false_positive, hash_functions, max_split_bins):
    """
    Return the correction factor for element split in multiple bins
    This ratio is how many time the bins have to grow to achieve the desired fpr
    """
    return math.log(1-(1-(1-false_positive)**max_split_bins)**(1/hash_functions))/math.log(1-false_positive**(1/hash_functions))


def estimate_params(cfg, simulated_bin_lens, groups_len):
    # Caculate filter sizes based on simulated bin lenghts
    params = {}
    for bin_length in simulated_bin_lens:

        # Estimate number of bins based on ranks and sequence sizes
        n_bins, max_split_bins = estimate_n_bins(bin_length, cfg.overlap_length, groups_len)
        n_bins = optimal_bins(n_bins)
        
        if n_bins <= 0:
            continue  # invalid bin_len

        max_elements = estimate_elements(bin_length, cfg.kmer_size, cfg.window_size)
        # Derive Bloom filter values from max_fp
        bf_params = derive_bf_params(max_elements,
                                     cfg.max_fp,
                                     0,
                                     cfg.hash_functions)

        if bf_params["hash_functions"] == 0:
            continue  # filter too small

        # Calculate correction ratio
        try:
            ratio = split_correction_ratio(bf_params["false_positive"], bf_params["hash_functions"], max_split_bins)
        except:
            # too many splits, ratio too high
            break

        params[bin_length] = bf_params
        params[bin_length]["n_bins"] = n_bins
        params[bin_length]["filter_size_bits"] = bf_params["size_bits"] * n_bins
        params[bin_length]["max_split_bins"] = max_split_bins
        params[bin_length]["corr_ratio"] = ratio
        params[bin_length]["corr_filter_size_bits"] = bf_params["size_bits"] * n_bins * ratio
 
    return params


def estimate_bin_length(cfg, seqinfo, tax):
    """
    Estimate best size of bin_length to split sequences. Always use --max-fp to estimate, even when --filter-size is provided
    """

    # number of simulations
    nsim = 2000

    # Generate dict with target groups and their total length to estimate params
    groups_len = {}
    if cfg.specialization:
        groups_len = seqinfo.seqinfo.groupby('specialization').sum().to_dict()['length']
    elif cfg.rank == "leaves":
        groups_len = seqinfo.seqinfo.groupby('taxid').sum().to_dict()['length']
    else:
        groups_len = pd.concat([seqinfo.seqinfo['taxid'].apply(lambda x: tax.get_rank(x, cfg.rank)), seqinfo.seqinfo['length']], axis=1).groupby('taxid').sum().to_dict()['length']

    # Biggest target group as max possible bin_len
    max_bin_len = max(groups_len.values())

    # Simulate possible ibf outcomes from biggest target group (max. bin len) by generating bin lenghts in rev. geometric space
    # Use geometric to have more simulations on bigger sizes
    simulated_bin_lens = map(round, np.geomspace(max_bin_len, 1, num=nsim, dtype=int))

    # simulated parameters
    # values are very different from true values with minimizers but proportions are a good measure
    params = estimate_params(cfg, simulated_bin_lens, groups_len)

    if not params:
        #could not estimate --bin-length, using max_bin_len
        return max_bin_len
    
    # Filter params for entires with 1.5x times the number of min. bins
    if cfg.faster:
        min_n_bins = min([v["n_bins"] for v in params.values()])
        params = dict(filter(lambda v: v[1]["n_bins"] <= optimal_bins(min_n_bins*1.5), params.items()))

    # Select best final corrected filter size
    # This value is just an approx. since it considers all possible kmers (minimizers will be way smaller and different)
    # When estimating with fixed filter_size, also use this as an trade-off metric for best parameters
    selected_bin_len = sorted(params, key=lambda k: params[k]["corr_filter_size_bits"])[0]

    #print('bin_length', 'hash_functions', 'n_bins', 'max_split_bins', 'corr_ratio', 'corr_filter_size_bits', sep="\t")
    #for p, v in params.items():
        #print(p, v)
        #print(p, v['hash_functions'], v['n_bins'], v['max_split_bins'], str("{0:.5f}".format(v['corr_ratio'])), str("{0:.2f}".format(bits2mb(v['corr_filter_size_bits']))), sep="\t")
    #print("")
    #print(params[sorted(params, key=lambda k: params[k]["corr_filter_size_bits"])[0]])

    return selected_bin_len


def run_taxsbp(seqinfo, bin_length, fragment_length, overlap_length, rank, specialization, ncbi_nodes_file, ncbi_merged_file, verbose, bins: Bins=None):
    taxsbp_params = {}

    taxsbp_params["input_table"] = seqinfo.seqinfo
    if bins is not None:
        taxsbp_params["update_table"] = bins.bins

    taxsbp_params["nodes_file"] = ncbi_nodes_file
    if ncbi_merged_file:
        taxsbp_params["merged_file"] = ncbi_merged_file

    taxsbp_params["bin_len"] = bin_length
    if fragment_length:
        taxsbp_params["fragment_len"] = fragment_length
        taxsbp_params["overlap_len"] = overlap_length

    if specialization:
        taxsbp_params["specialization"] = specialization
        taxsbp_params["bin_exclusive"] = specialization
    else:  # either species,genus ... or "leaves"
        taxsbp_params["bin_exclusive"] = rank

    taxsbp_params["silent"] = False if verbose else True

    return Bins(taxsbp_ret=taxsbp.taxsbp.pack(**taxsbp_params), use_specialization=True if specialization else False)


def derive_bf_params(elements, false_positive, size_bits, hash_functions):
    """
    given a number of elements and (false_positive or size_bits) derive minimal values
    returns a dict with  {"size_bits", "false_positive", "hash_functions"}
    """
    max_hash_functions = 5

    def ratio_from_size_elements(size, elements):
        return size / elements

    def ratio_from_hashf_fp(hashf, fp):
        return -hashf / math.log(1 - math.exp(math.log(fp) / hashf))

    def hashf_from_ratio(ratio):
        n = round(math.log(2) * ratio)
        return n if n <= max_hash_functions else max_hash_functions

    def hashf_from_fp(fp):
        n = -round(math.log2(fp))
        return n if n <= max_hash_functions else max_hash_functions

    if size_bits:
        final_size_bits = size_bits
        r = ratio_from_size_elements(final_size_bits, elements)
        final_hash_functions = hashf_from_ratio(r) if not hash_functions else hash_functions
        final_false_positive = 1 if not final_hash_functions else math.pow(1 - math.exp(-final_hash_functions / r), final_hash_functions)
    elif false_positive:
        final_false_positive = false_positive
        final_hash_functions = hashf_from_fp(final_false_positive) if not hash_functions else hash_functions
        r = 0 if not final_hash_functions else ratio_from_hashf_fp(final_hash_functions, final_false_positive)
        final_size_bits = math.ceil(elements * r)
    else:
        final_size_bits = size_bits
        final_false_positive = false_positive
        final_hash_functions = hash_functions

    return {"size_bits": final_size_bits,
            "false_positive": final_false_positive,
            "hash_functions": final_hash_functions}
