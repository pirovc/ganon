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
        print_log("ERROR: No valid input files found")
        return False

    # Set db prefixes
    db_prefix = {prefix:cfg.db_prefix + "." + prefix for prefix in  ["ibf","map","tax","gnn"]}  
    
    # Set temporary working folder 
    tmp_output_folder = cfg.db_prefix + "_tmp/"
    if not set_tmp_folder(tmp_output_folder): return False

    # Set up taxonomy
    ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(cfg.taxdump_file, tmp_output_folder, cfg.quiet)
    
    tx = time.time()
    print_log("Parsing taxonomy", cfg.quiet)
    tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    seqinfo = SeqInfo()
    if cfg.seq_info_file:
        tx = time.time()
        print_log("Parsing seq-info-file", cfg.quiet)
        seqinfo.parse_seq_info_file(cfg.seq_info_file, parse_specialization=True if cfg.specialization=="custom" else False)
        
        # check if file has specialization col wiht valid values
        if cfg.specialization=="custom" and seqinfo.any_null("specialization"):
            print_log(" - Skipping custom specialization. Invalid values in 4th column of the --seq-info-file " + cfg.seq_info_file, cfg.quiet)
            cfg.specialization=""

        print_log(" - "  + str(seqinfo.size()) + " sequence entries in the --seq-info-file " + cfg.seq_info_file, cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        tx = time.time()
        print_log("Extracting sequence identifiers", cfg.quiet)
        parse_seqids(seqinfo, input_files, cfg.specialization, get_length=False)
        print_log(" - "  + str(seqinfo.size()) + " sequence headers successfully retrieved from " + str(len(input_files)) + " input file(s)" , cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
        retrieve_seqinfo(seqinfo, tmp_output_folder, input_files, cfg)
        if cfg.write_seq_info_file: seqinfo.write(cfg.db_prefix+".seqinfo.txt")

    # Convert cols data types
    seqinfo.convert()

    # check sequences compared to bins
    added_seqids, _, _ = check_updated_seqids(set(seqinfo.get_seqids()), set())
    # Ignore removed sequences if not doing complete update
    print_log("Build: adding " + str(len(added_seqids)) + " sequences", cfg.quiet)
    print_log("", cfg.quiet)

    # Set bin length
    if cfg.bin_length: # user defined
        bin_length = cfg.bin_length
    else:
        tx = time.time()
        print_log("Calculating best bin length", cfg.quiet)
        bin_length, approx_size, n_bins = estimate_bin_len_size(cfg, seqinfo, tax)
        if bin_length<=0: 
            bin_length=1000000
            print_log("WARNING: could not estimate bin length, using default of " + str(bin_length) + "bp")
        else:
            print_log(" - bin length: " + str(bin_length) + "bp (approx: " + str(n_bins) + " bins / " + str("{0:.2f}".format(approx_size)) + "MB)", cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # Set fragment length
    if cfg.fragment_length==-1: # if ==-1 set default
        fragment_length = bin_length - cfg.overlap_length
    elif cfg.fragment_length==0: # if ==0 deactivate
        fragment_length = 0
    else: # user input
        fragment_length = cfg.fragment_length - cfg.overlap_length 

    tx = time.time()
    print_log("Running taxonomic clustering (TaxSBP)", cfg.quiet)
    bins = run_taxsbp(seqinfo, bin_length, fragment_length, cfg.overlap_length, cfg.rank, cfg.specialization, ncbi_nodes_file, ncbi_merged_file)
    # bin statistics
    actual_number_of_bins = bins.get_number_of_bins()
    optimal_number_of_bins = optimal_bins(actual_number_of_bins)
    max_length_bin = bins.get_max_bin_length()
    max_kmer_count = max_length_bin-cfg.kmer_size+1 # aproximate number of unique k-mers by just considering that they are all unique
    print_log(" - " + str(actual_number_of_bins) + " bins created", cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

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
    if cfg.specialization: tax.add_nodes(bins.get_specialization_taxid(), cfg.specialization) 
    tax.write(db_prefix["tax"])

    if cfg.specialization and cfg.rank!="leaves":
        print_log(" - --rank is set to leaves when using specialization values", cfg.quiet)
        cfg.rank="leaves"
    # Write .gnn file
    print_log(" - " + db_prefix["gnn"], cfg.quiet)
    gnn = Gnn(kmer_size=cfg.kmer_size, 
            hash_functions=cfg.hash_functions, 
            number_of_bins=actual_number_of_bins, 
            rank=cfg.rank,
            bin_length=bin_length,
            fragment_length=fragment_length,
            overlap_length=cfg.overlap_length,
            bins=bins.get_list())
    gnn.write(db_prefix["gnn"])
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    print_log("Building index (ganon-build)", cfg.quiet)
    # define bloom filter size based on given false positive
    MBinBits = 8388608
    print_log(" - max unique " + str(cfg.kmer_size) +  "-mers: " + str(max_kmer_count), cfg.quiet)
    if not cfg.fixed_bloom_size:
        bin_size_bits = math.ceil(-(1/((1-cfg.max_fp**(1/float(cfg.hash_functions)))**(1/float(cfg.hash_functions*max_kmer_count))-1)))   
        print_log(" - IBF calculated size with fp<=" + str(cfg.max_fp) + ": " + str("{0:.2f}".format((bin_size_bits*optimal_number_of_bins)/MBinBits)) + "MB (" + str(bin_size_bits) + " bits/bin * " + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])", cfg.quiet)
    else:
        bin_size_bits = math.ceil((cfg.fixed_bloom_size * MBinBits)/optimal_number_of_bins);
        estimated_max_fp = (1-((1-(1/float(bin_size_bits)))**(cfg.hash_functions*max_kmer_count)))**cfg.hash_functions
        print_log(" - IBF calculated max. fp with size=" + str(cfg.fixed_bloom_size) + "MB: " + str("{0:.2f}".format(estimated_max_fp) + " ("  + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])"), cfg.quiet)

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
                                    "--filter-size-bits " + str(bin_size_bits*optimal_number_of_bins) if cfg.max_fp else "--filter-size " + str(cfg.fixed_bloom_size),
                                    "--kmer-size " + str(cfg.kmer_size),
                                    "--hash-functions " + str(cfg.hash_functions),
                                    "--threads " + str(cfg.threads),
                                    "--output-filter-file " + db_prefix["ibf"],
                                    "--verbose" if cfg.verbose else "",
                                    "--quiet" if cfg.quiet else "",
                                    "--n-refs " + str(cfg.n_refs) if cfg.n_refs is not None else "",
                                    "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                    "--reference-files " + ",".join([file for file in input_files]) if input_files else "",
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
        print_log("ERROR: No valid input files found")
        return False

    # Set db prefixes
    db_prefix = {prefix:cfg.db_prefix + "." + prefix for prefix in  ["ibf","map","tax","gnn"]}  
    
    # Set temporary working folder (current or new output)
    tmp_output_folder = cfg.output_db_prefix + "_tmp/" if cfg.output_db_prefix else cfg.db_prefix + "_tmp/"
    if not set_tmp_folder(tmp_output_folder): return False

    # Load .gnn file   
    gnn = Gnn(file=db_prefix["gnn"])

    # load bins
    bins = Bins(taxsbp_ret=gnn.bins)

    seqinfo = SeqInfo()
    if cfg.seq_info_file:
        tx = time.time()
        print_log("Parsing --seq-info-file", cfg.quiet)
        seqinfo.parse_seq_info_file(cfg.seq_info_file, parse_specialization=True if cfg.specialization=="custom" else False)
        # check if file has specialization col wiht valid values
        if cfg.specialization=="custom" and seqinfo.any_null("specialization"):
            print_log(" - Skipping custom specialization. Invalid values in 4th column of the --seq-info-file " + cfg.seq_info_file, cfg.quiet)
            cfg.specialization=""
        print_log(" - "  + str(seqinfo.size()) + " sequence entries in the --seq-info-file " + cfg.seq_info_file, cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        tx = time.time()
        print_log("Extracting sequence identifiers", cfg.quiet)
        parse_seqids(seqinfo, input_files, cfg.specialization, get_length=False)
        print_log(" - "  + str(seqinfo.size()) + " sequence headers successfully retrieved from " + str(len(input_files)) + " input file(s)" , cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    # check sequences compared to bins
    added_seqids, removed_seqids, kept_seqids = check_updated_seqids(set(seqinfo.get_seqids()), set(bins.get_seqids()))
    # Ignore removed sequences if not doing complete update
    if cfg.update_complete: 
        print_log("Update: adding " + str(len(added_seqids)) + " sequences, removing " + str(len(removed_seqids)) + " sequences, keeping " + str(len(kept_seqids)) + " sequences", cfg.quiet)
    else:
        removed_seqids=[]
        print_log("Update: adding " + str(len(added_seqids)) + " sequences, ignoring " + str(len(kept_seqids)) + " repeated sequences", cfg.quiet)
    print_log("", cfg.quiet)

    if not added_seqids and not removed_seqids:
        print_log("Nothing to update", cfg.quiet)
        return False

    if cfg.update_complete:           
        # Remove already included seqids to just retrieve information for added sequences
        seqinfo.remove_seqids(kept_seqids | removed_seqids)
    else:
        # Remove seqids already present in the current version (repeated entries)
        seqinfo.remove_seqids(kept_seqids)

    # load seqinfo file with data (after removing ids)
    if not cfg.seq_info_file: 
        retrieve_seqinfo(seqinfo, tmp_output_folder, input_files, cfg)
        if cfg.write_seq_info_file: seqinfo.write(cfg.output_db_prefix+".seqinfo.txt")

    # Convert cols data types
    seqinfo.convert()
    
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
    # Compability to nomenclature before 0.3.5
    if gnn.rank=="taxid" or gnn.rank=="assembly": gnn.rank="leaves"  
    updated_bins = run_taxsbp(seqinfo, gnn.bin_length, gnn.fragment_length, gnn.overlap_length, gnn.rank, cfg.specialization, ncbi_nodes_file, ncbi_merged_file, bins=bins)
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
    if cfg.specialization: tax.add_nodes(updated_bins.get_specialization_taxid(), cfg.specialization) 
    
    # Load old .tax file into new taxonomy
    tax.merge(Tax([db_prefix["tax"]]))
    # Write .tax file
    tax.write(cfg.output_db_prefix + ".tax" if cfg.output_db_prefix else db_prefix["tax"])
    # TODO - remove entries from .tax from removed entries of the db

    # merge updated and old bins together
    bins.merge(updated_bins)
    
    # Write .gnn file
    print_log(" - " + cfg.output_db_prefix + ".gnn" if cfg.output_db_prefix else db_prefix["gnn"], cfg.quiet)
    gnn.bins = bins.get_list() # save updated bins
    gnn.number_of_bins=bins.get_number_of_bins() # add new bins count
    gnn.write(cfg.output_db_prefix + ".gnn" if cfg.output_db_prefix else db_prefix["gnn"])

    # Recreate .map file based on the new bins
    print_log(" - " + cfg.output_db_prefix + ".map" if cfg.output_db_prefix else db_prefix["map"], cfg.quiet)
    bins.write_map_file(cfg.output_db_prefix + ".map" if cfg.output_db_prefix else db_prefix["map"], use_specialization=True if cfg.specialization else False)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    tx = time.time()
    print_log("Updating index (ganon-build)", cfg.quiet)

    # Write aux. file for ganon
    # This file has to contain all new sequences
    # in case of update_complete, 
    acc_bin_file = tmp_output_folder + "acc_bin.txt"
    
    if cfg.update_complete:
        # all sequences from the bins with added/removed sequences should be written
        bins.write_acc_bin_file(acc_bin_file, new_binids | updated_binids)
        # If all sequences of a bin were removed and no new sequence added
        # insert a dummy entry for ganon-build to clear the bin
        if removed_binids: 
            with open(acc_bin_file,"a") as abf:
                for b in removed_binids:
                    print(0,0,0,b,sep="\t",file=abf)

    else:
        # Only new sequences (updated_bins) either on old or new binids
        updated_bins.write_acc_bin_file(acc_bin_file)

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
                                    "--seqid-bin-file " + acc_bin_file,
                                    "--output-filter-file " + tmp_db_prefix_ibf,
                                    "--threads " + str(cfg.threads),                                
                                    "--verbose" if cfg.verbose else "",
                                    "--quiet" if cfg.quiet else "",
                                    "--n-refs " + str(cfg.n_refs) if cfg.n_refs is not None else "",
                                    "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                    "--reference-files " + ",".join([file for file in input_files]) if input_files else "",
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

def parse_seqids(seqinfo, input_files, specialization, get_length: bool):
    for file in input_files:
        if get_length:
            # cat | zcat | gawk -> compability with osx
            run_get = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{if (seqlen){{print seqlen}}; printf substr($1,2)\"\\t\";seqlen=0;next;}} {{seqlen+=length($0)}}END{{print seqlen}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get, print_stderr=False, shell=True)
            parsed_stdout = pd.read_csv(StringIO(stdout), sep="\t", header=None, names=['seqid', 'length'])
        else:
            # cat | zcat | gawk -> compability with osx
            run_get = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{print substr($1,2)}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get, print_stderr=False, shell=True)
            parsed_stdout = pd.read_csv(StringIO(stdout), header=None, names=['seqid'])
        if specialization=="file":
            parsed_stdout["specialization"] = os.path.basename(file)
        elif specialization=="sequence":
            parsed_stdout["specialization"] = parsed_stdout["seqid"]
        seqinfo.append(parsed_stdout)

    if get_length:
        # Drop rows with zero length
        seqinfo.drop_zeros(col='length')

def parse_eutils(seqinfo, tmp_output_folder, path_exec_get_seq_info, skip_len_taxid=False, get_assembly=False):
    seqid_file = tmp_output_folder + "seqids.txt"
    seqinfo.write_seqid_file(seqid_file)
    run_get_seq_info_cmd = '{0} -k -r -i {1} {2} {3}'.format(
                                path_exec_get_seq_info,
                                seqid_file,
                                "-a" if get_assembly else "",
                                "-s" if skip_len_taxid else "")
    stdout, stderr = run(run_get_seq_info_cmd, print_stderr=True, exit_on_error=False)
    
    # always return all entries in the same order (-k)
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

    # drop failed entries for lenght or taxid, assembly is NaN if no specialization
    seqinfo.dropna(subset=["length","taxid"])

def retrieve_seqinfo(seqinfo, tmp_output_folder, input_files, cfg):
    # Max. # of sequences to use eutils as auto mode
    max_seqs_eutils = 50000
    # initialize total count
    seqid_total_count = seqinfo.size()

    # Define method to use
    if "auto" in cfg.seq_info_mode:
        if seqid_total_count>max_seqs_eutils: 
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
        parse_eutils(seqinfo, tmp_output_folder, cfg.path_exec['get_seq_info'], skip_len_taxid=False, get_assembly=True if cfg.specialization=="assembly" else False)
        print_log(" - " + str(seqinfo.size()) + " sequences successfully retrieved", cfg.quiet)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        # Retrieve seq. ids and lengths
        tx = time.time()
        print_log("Extracting sequence lengths", cfg.quiet)
        seqinfo.clear() # Clear seqinfo (get seqids again with length)
        parse_seqids(seqinfo, input_files, cfg.specialization, get_length=True)
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
        # Check if retrieved taxids are the same as number of inputs, reset counter
        if seqinfo.size() < seqid_total_count:
            print_log(" - could not retrieve taxid for " + str(seqid_total_count - seqinfo.size()) + " accessions", cfg.quiet)
            seqid_total_count = seqinfo.size()
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

        if cfg.specialization=="assembly":
            tx = time.time()
            print_log("Retrieving assembly information from NCBI E-utils", cfg.quiet)
            parse_eutils(seqinfo, tmp_output_folder, cfg.path_exec['get_seq_info'], skip_len_taxid=True, get_assembly=True)
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
        tmp_seqid_taxids = pd.read_csv(acc2txid, sep='\t', header=None, skiprows=1, usecols=[1,2], names=['seqid','taxid'], converters={'seqid':lambda x: x if x in unique_seqids else ""}, dtype={'taxid': 'str'})
        tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['seqid']!=""] #keep only seqids used
        tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['taxid']!="0"] # filter out taxid==0
        # save count to return
        count_acc2txid[acc2txid] = tmp_seqid_taxids.shape[0]    
        # merge taxid retrieved based on seqid
        seqinfo.join(tmp_seqid_taxids, "taxid")
        del tmp_seqid_taxids
        #if already found all seqids no need to parse all files till the end)
        if sum(count_acc2txid.values()) == len(unique_seqids): 
            break 

    return count_acc2txid

def ibf_size_mb(bin_len, n_bins, max_fp, hash_functions, kmer_size):
    return (math.ceil(-(1/((1-max_fp**(1/float(hash_functions)))**(1/float(hash_functions*(bin_len-kmer_size+1)))-1)))*optimal_bins(n_bins))/8388608

def approx_n_bins(bin_len, overlap_len, groups_len): 
    frag_len=bin_len-overlap_len
    n_bins = sum([math.ceil(math.ceil(l/(frag_len-overlap_len))/(bin_len/(frag_len+overlap_len))) for l in groups_len.values()])
    return n_bins

def optimal_bins(n): 
    #return optimal number of bins for the IBF (multiples of 64)
    return (math.floor(n/64)+1)*64 

def estimate_bin_len_size(cfg, seqinfo, tax):
    # Simulate bins that will be created by taxsbp with many bin lenghts
    # Select the best trade-off on size and n. of bins

    # Generate dict with groups:total_length
    groups_len = {}
    if cfg.specialization:
        groups_len = seqinfo.seqinfo.groupby('specialization').sum().to_dict()['length']
    elif cfg.rank=="leaves":
        groups_len = seqinfo.seqinfo.groupby('taxid').sum().to_dict()['length']
    else:
        groups_len = pd.concat([seqinfo.seqinfo['taxid'].apply(lambda x: tax.get_rank(x, cfg.rank)), seqinfo.seqinfo['length']], axis=1).groupby('taxid').sum().to_dict()['length']

    # Set limits
    # fixed start size (too low generates few possible cases for analysis)
    min_bin_len = 500
    # Biggest group as max
    max_bin_len = max(groups_len.values())
    # Try to find min. size by simulating points in geometric space 
    # between min. and max. bin length. Use geometric to have more simulations on smaller sizes
    # Necessary to calculate instead of getting first lowest
    # since it may be a local minimum
    bin_lens = np.geomspace(min_bin_len, max_bin_len, num=300)

    # Caculate filter sizes based on bin_lens
    filter_sizes = np.array([ibf_size_mb(b, approx_n_bins(b, cfg.overlap_length, groups_len), cfg.max_fp, cfg.hash_functions, cfg.kmer_size) for b in bin_lens])
    # keep only valid positive entries to define min. filter size
    idx_above_min = filter_sizes>0
    filter_sizes = filter_sizes[idx_above_min]
    bin_lens = bin_lens[idx_above_min]
    min_filter_size=filter_sizes.min()

    print_log(" - Approx. min. size possible: " + str("{0:.2f}".format(min_filter_size)) + "MB", cfg.quiet)
    # Define max size
    # if none defined or too small, Define the max as 1.5 time size of the min
    if cfg.max_bloom_size is None:
        max_filter_size = min_filter_size*1.5
    elif cfg.max_bloom_size<min_filter_size:
        max_filter_size = min_filter_size*1.5
        print_log(" - --max-bloom-size " + str(cfg.max_bloom_size) + "MB is too small, using max. default (1.5x min.): " + str("{0:.2f}".format(max_filter_size)) + "MB", cfg.quiet)
    else:
        max_filter_size = cfg.max_bloom_size
    
    # keep only valid points below max_filter_size
    idx_below_max = filter_sizes<=max_filter_size

    # If more than one valid point
    if sum(idx_below_max)>1: 
        # reduce space in between min. and max. filter size to get better
        bin_lens = np.linspace(bin_lens[idx_below_max].min(), bin_lens[idx_below_max].max(), num=300)
        # Estimate n_bins
        n_bins = [approx_n_bins(b, cfg.overlap_length, groups_len) for b in bin_lens]
        filter_sizes = [ibf_size_mb(b, n_bins[i], cfg.max_fp, cfg.hash_functions, cfg.kmer_size) for i,b in enumerate(bin_lens)]
        # Get value with min. number of bins from this distribution
        idx_min = np.where(n_bins == np.amin(n_bins))[0][0]
        return int(bin_lens[idx_min]), filter_sizes[idx_min], n_bins[idx_min]    
    else:
        return 0,0,0

def run_taxsbp(seqinfo, bin_length, fragment_length, overlap_length, rank, specialization, ncbi_nodes_file, ncbi_merged_file, bins: Bins=None):
    taxsbp_params={}

    taxsbp_params["input_table"] = seqinfo.to_csv()
    if bins is not None:
        taxsbp_params["update_table"] = bins.to_csv()

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
    else: # either species,genus ... or "leaves"
        taxsbp_params["bin_exclusive"] =  rank

    return Bins(taxsbp_ret=taxsbp.taxsbp.pack(**taxsbp_params))