import time, math
import pandas as pd
import numpy as np
import taxsbp.taxsbp
from src.ganon.bins import Bins
from src.ganon.gnn import Gnn
from src.ganon.seqinfo import SeqInfo
from src.ganon.tax import Tax
from src.ganon.util import *

def build(cfg):
    # validate input files
    input_files, input_files_from_directory = validate_input_files(cfg)
    if len(input_files)==0 and len(input_files_from_directory)==0:
        print_log("No valid input files found")
        sys.exit(1)

    # Set db prefixes
    db_prefix = {prefix:cfg.db_prefix + "." + prefix for prefix in  ["ibf","map","tax","gnn"]}  
    
    # Set temporary working folder 
    tmp_output_folder = cfg.db_prefix + "_tmp/"
    set_tmp_folder(tmp_output_folder)

    # Set assembly mode
    use_assembly=True if cfg.rank=="assembly" else False

    # Set up taxonomy
    ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(cfg, tmp_output_folder)
    
    tx = time.time()
    print_log("Parsing taxonomy")
    tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

    # Load seqids and generate seqinfo
    if cfg.seq_info_file:
        seqinfo = load_seqids(seq_info_file=cfg.seq_info_file)
    else:
        seqinfo = load_seqids(files=input_files + input_files_from_directory) 
        load_seqinfo(tmp_output_folder, seqinfo, cfg.path_exec, cfg.seq_info_mode, use_assembly)

    # Set bin length
    if cfg.bin_length: # user defined
        bin_length = cfg.bin_length
    else:
        tx = time.time()
        print_log("Calculating best bin length")
        bin_length, approx_size, n_bins = estimate_bin_len_size(cfg, seqinfo, tax, use_assembly)
        if bin_length<=0: 
            bin_length=1000000
            print_log(" - could not estimate bin length, using default of " + str(bin_length) + "bp")
        else:
            print_log(" - bin length: " + str(bin_length) + "bp (approx: " + str(n_bins) + " bins / " + str("{0:.2f}".format(approx_size)) + "MB)")
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

    # Set fragment length
    if cfg.fragment_length==-1: # if ==-1 set default
        fragment_length = bin_length - cfg.overlap_length
    elif cfg.fragment_length==0: # if ==0 deactivate
        fragment_length = 0
    else: # user input
        fragment_length = cfg.fragment_length - cfg.overlap_length 

    tx = time.time()
    print_log("Running taxonomic clustering (TaxSBP)")
    taxsbp_params={}
    taxsbp_params["nodes_file"] = ncbi_nodes_file
    taxsbp_params["bin_len"] = bin_length
    if use_assembly:
        taxsbp_params["bin_exclusive"] = "assembly"
    elif cfg.rank=="taxid":
        taxsbp_params["bin_exclusive"] = "leaves"
    else:
        taxsbp_params["bin_exclusive"] =  cfg.rank
    if ncbi_merged_file: taxsbp_params["merged_file"] = ncbi_merged_file
    if fragment_length: 
        taxsbp_params["fragment_len"] = fragment_length
        taxsbp_params["overlap_len"] = cfg.overlap_length
    if use_assembly: taxsbp_params["specialization"] = "assembly"
    taxsbp_params["input_table"] = seqinfo.get_csv()
    bins = Bins(taxsbp_ret=taxsbp.taxsbp.pack(**taxsbp_params))
    del taxsbp_params
    # bin statistics
    actual_number_of_bins = bins.get_number_of_bins()
    optimal_number_of_bins = optimal_bins(actual_number_of_bins)
    max_length_bin = bins.get_max_bin_length()
    max_kmer_count = max_length_bin-cfg.kmer_size+1 # aproximate number of unique k-mers by just considering that they are all unique
    print_log(" - " + str(actual_number_of_bins) + " bins created")
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

    tx = time.time()
    print_log("Building database files")
    
    # Write .map file
    print_log(" - " + db_prefix["map"])
    bins.write_map_file(db_prefix["map"], use_assembly)

    # Write .tax file
    print_log(" - " + db_prefix["tax"])
    tax.filter(bins.get_taxids()) # filter only used taxids
    if use_assembly: tax.add_nodes(bins.get_specialization_taxid(), "assembly") # add assembly nodes
    tax.write(db_prefix["tax"])

    # Write .gnn file
    print_log(" - " + db_prefix["gnn"])
    gnn = Gnn(kmer_size=cfg.kmer_size, 
            hash_functions=cfg.hash_functions, 
            number_of_bins=actual_number_of_bins, 
            rank=cfg.rank,
            bin_length=bin_length,
            fragment_length=fragment_length,
            overlap_length=cfg.overlap_length,
            bins=bins.get_list())
    gnn.write(db_prefix["gnn"])
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

    print_log("Building index (ganon-build)")
    # define bloom filter size based on given false positive
    MBinBits = 8388608
    print_log(" - max unique " + str(cfg.kmer_size) +  "-mers: " + str(max_kmer_count))
    if not cfg.fixed_bloom_size:
        bin_size_bits = math.ceil(-(1/((1-cfg.max_fp**(1/float(cfg.hash_functions)))**(1/float(cfg.hash_functions*max_kmer_count))-1)))   
        print_log(" - IBF calculated size with fp<=" + str(cfg.max_fp) + ": " + str("{0:.2f}".format((bin_size_bits*optimal_number_of_bins)/MBinBits)) + "MB (" + str(bin_size_bits) + " bits/bin * " + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])")
    else:
        bin_size_bits = math.ceil((cfg.fixed_bloom_size * MBinBits)/optimal_number_of_bins);
        estimated_max_fp = (1-((1-(1/float(bin_size_bits)))**(cfg.hash_functions*max_kmer_count)))**cfg.hash_functions
        print_log(" - IBF calculated max. fp with size=" + str(cfg.fixed_bloom_size) + "MB: " + str("{0:.2f}".format(estimated_max_fp) + " ("  + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])"))

    # Write aux. file for ganon
    acc_bin_file = tmp_output_folder + "acc_bin.txt"
    bins.write_acc_bin_file(acc_bin_file)

    run_ganon_build_cmd = " ".join([cfg.path_exec['build'],
                                    "--seqid-bin-file " + acc_bin_file,
                                    "--filter-size-bits " + str(bin_size_bits*optimal_number_of_bins) if cfg.max_fp else "--filter-size " + str(cfg.fixed_bloom_size),
                                    "--kmer-size " + str(cfg.kmer_size),
                                    "--hash-functions " + str(cfg.hash_functions),
                                    "--threads " + str(cfg.threads),
                                    "--output-filter-file " + db_prefix["ibf"],
                                    "--verbose" if cfg.verbose else "",
                                    "--n-refs " + str(cfg.n_refs) if cfg.n_refs is not None else "",
                                    "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                    "--reference-files " + ",".join([file for file in input_files]) if input_files else "",
                                    "--directory-reference-files " + cfg.input_directory if cfg.input_directory else "",
                                    "--extension " + cfg.input_extension if cfg.input_extension else ""])
    stdout, stderr = run(run_ganon_build_cmd, print_stderr=True)

    # Delete temp files
    rm_tmp_folder(tmp_output_folder)

def update(cfg):
    tx = time.time()
    # validate input files
    input_files, input_files_from_directory = validate_input_files(cfg)
    if len(input_files)==0 and len(input_files_from_directory)==0:
        print_log("No valid input files found")
        sys.exit(1)

    # Set db prefixes
    db_prefix = {prefix:cfg.db_prefix + "." + prefix for prefix in  ["ibf","map","tax","gnn"]}  
    
    # Set temporary working folder (current or new output)
    tmp_output_folder = cfg.output_db_prefix + "_tmp/" if cfg.output_db_prefix else cfg.db_prefix + "_tmp/"
    set_tmp_folder(tmp_output_folder)

    # Load .gnn file   
    gnn = Gnn(file=db_prefix["gnn"])
    # Set assembly mode
    use_assembly=True if gnn.rank=="assembly" else False

    # load bins
    bins = Bins(taxsbp_ret=gnn.bins)

    # Load seqids and generate seqinfo
    if cfg.seq_info_file:
        seqinfo = load_seqids(seq_info_file=cfg.seq_info_file)
    else:
        seqinfo = load_seqids(files=input_files + input_files_from_directory) 

    # check sequences compared to bins
    added_seqids, removed_seqids, kept_seqids = check_updated_seqids(set(seqinfo.get_seqids()), set(bins.get_seqids()))
    # Ignore removed sequences if not doing complete update
    if cfg.update_complete: 
        print_log("Update: adding " + str(len(added_seqids)) + " sequences, removing " + str(len(removed_seqids)) + " sequences, keeping " + str(len(kept_seqids)) + " sequences")
    else:
        removed_seqids=[]
        print_log("Update: adding " + str(len(added_seqids)) + " sequences, ignoring " + str(len(kept_seqids)) + " repeated sequences")
    print_log("")

    if not added_seqids and not removed_seqids:
        print_log("Nothing to update.")
        sys.exit(0)

    if cfg.update_complete:           
        # Remove already included seqids to just retrieve information for added sequences
        seqinfo.remove_seqids(kept_seqids | removed_seqids)
    else:
        # Remove seqids already present in the current version (repeated entries)
        seqinfo.remove_seqids(kept_seqids)

    # load seqinfo file with data (after removing ids)
    if not cfg.seq_info_file: load_seqinfo(tmp_output_folder, seqinfo, cfg.path_exec, cfg.seq_info_mode, use_assembly)

    # save set of current binids
    previous_binids = set(bins.get_binids())
    # remove seqids from bins if performing update complete
    if cfg.update_complete and removed_seqids:
        bins.remove_seqids(removed_seqids)
    # save set of kept binids after removal
    kept_binids = set(bins.get_binids())

    # Set up taxonomy files
    ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(cfg, tmp_output_folder)

    tx = time.time()
    print_log("Running taxonomic clustering (TaxSBP)")
    taxsbp_params={}
    taxsbp_params["update_table"] = bins.get_csv()
    taxsbp_params["nodes_file"] = ncbi_nodes_file
    taxsbp_params["bin_len"] = gnn.bin_length
    if use_assembly:
        taxsbp_params["bin_exclusive"] = "assembly"
    elif gnn.rank=="taxid":
        taxsbp_params["bin_exclusive"] = "leaves"
    else:
        taxsbp_params["bin_exclusive"] =  gnn.rank
    if ncbi_merged_file: taxsbp_params["merged_file"] = ncbi_merged_file
    if gnn.fragment_length: 
        taxsbp_params["fragment_len"] = gnn.fragment_length
        taxsbp_params["overlap_len"] = gnn.overlap_length
    if use_assembly: taxsbp_params["specialization"] = "assembly"
    taxsbp_params["input_table"] = seqinfo.get_csv()
    updated_bins = Bins(taxsbp_ret=taxsbp.taxsbp.pack(**taxsbp_params))
    # bin statistics
    taxsbp_binids = set(updated_bins.get_binids())
    removed_binids = previous_binids.difference(kept_binids | taxsbp_binids)
    new_binids = taxsbp_binids.difference(previous_binids)
    updated_binids = taxsbp_binids.intersection(previous_binids)
    print_log(" - " + str(len(new_binids)) + " bins added, " + str(len(updated_binids)) + " bins updated, " + str(len(removed_binids)) + " bins removed")
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

    tx = time.time()
    print_log("Updating database files")
    # load new taxonomy
    print_log(" - " + cfg.output_db_prefix + ".tax" if cfg.output_db_prefix else db_prefix["tax"])
    tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
    # Update and write .tax file
    tax.filter(updated_bins.get_taxids()) # filter only used taxids
    if use_assembly: tax.add_nodes(updated_bins.get_specialization_taxid(), "assembly") # add assembly nodes
    # Load old .tax file into new taxonomy
    tax.merge(Tax([db_prefix["tax"]]))
    # Write .tax file
    tax.write(cfg.output_db_prefix + ".tax" if cfg.output_db_prefix else db_prefix["tax"])
    # TODO - remove entries from .tax from removed entries of the db

    # merge updated and old bins together
    bins.merge(updated_bins)
    
    # Write .gnn file
    print_log(" - " + cfg.output_db_prefix + ".gnn" if cfg.output_db_prefix else db_prefix["gnn"])
    gnn.bins = bins.get_list() # save updated bins
    gnn.number_of_bins=bins.get_number_of_bins() # add new bins count
    gnn.write(cfg.output_db_prefix + ".gnn" if cfg.output_db_prefix else db_prefix["gnn"])

    # Recreate .map file based on the new bins
    print_log(" - " + cfg.output_db_prefix + ".map" if cfg.output_db_prefix else db_prefix["map"])
    bins.write_map_file(cfg.output_db_prefix + ".map" if cfg.output_db_prefix else db_prefix["map"], use_assembly)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

    tx = time.time()
    print_log("Updating index (ganon-build)")

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

    # Temporary output filter 
    tmp_db_prefix_ibf = tmp_output_folder + "ganon.ibf"
    run_ganon_build_cmd = " ".join([cfg.path_exec['build'],
                                    "--update-filter-file " + db_prefix["ibf"],
                                    "--seqid-bin-file " + acc_bin_file,
                                    "--output-filter-file " + tmp_db_prefix_ibf,
                                    "--threads " + str(cfg.threads),                                
                                    "--verbose" if cfg.verbose else "",
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
    
def check_updated_seqids(new_seqids, old_seqids):
    # remove repeated from old bins
    added_seqids = new_seqids.difference(old_seqids)
    removed_seqids = old_seqids.difference(new_seqids)
    kept_seqids = old_seqids.difference(removed_seqids)

    return added_seqids, removed_seqids, kept_seqids


def load_seqids(files: list=[], seq_info_file: str=None):
    tx = time.time()
    print_log("Extracting sequence identifiers")
    # Load or create seqinfo
    if seq_info_file is not None: # file already provided 
        seqinfo = SeqInfo(seq_info_file=seq_info_file)
    else: # retrieve info
        # Count number of input sequences to define method or retrieve accessions for forced eutils
        seqinfo = SeqInfo()
        seqinfo.parse_seqid(files)
    print_log(" - "  + str(seqinfo.size()) + " sequences successfully retrieved")
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")
    return seqinfo

def load_seqinfo(tmp_output_folder, seqinfo, path_exec, seq_info_mode, use_assembly):
    
    # Max. # of sequences to use eutils as auto mode
    max_seqs_eutils = 50000
    # default accession2taxid files
    default_acc2txid = ["nucl_gb", "nucl_wgs"]
    # initialize total count
    seqid_total_count = seqinfo.size()

    # Define method to use
    if seq_info_mode[0]=="auto" and seqid_total_count>max_seqs_eutils: 
        seq_info_mode = default_acc2txid
    else:
        seq_info_mode = ["eutils"]

    if seq_info_mode[0]=="eutils":
        tx = time.time()
        print_log("Retrieving sequence information from NCBI E-utils")
        seqid_file = tmp_output_folder + "seqids.txt"
        seqinfo.write_seqid_file(seqid_file)
        seqinfo.parse_ncbi_eutils(seqid_file, path_exec['get_len_taxid'], skip_len_taxid=False, get_assembly=True if use_assembly else False)
        print_log(" - " + str(seqinfo.size()) + " sequences successfully retrieved")
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")
    else:
        # acc2taxid - offline mode
        acc2txid_options = ["nucl_gb","nucl_wgs","nucl_est","nucl_gss","pdb","prot","dead_nucl","dead_wgs","dead_prot"]

        # Retrieve seq. lengths
        tx = time.time()
        print_log("Extracting sequence lengths")
        seqinfo.parse_seqid_length(input_files)
        print_log(" - " + str(seqinfo.size()) + " sequences successfully retrieved")
        # Check if retrieved lengths are the same as number of inputs, reset counter
        if seqinfo.size() < seqid_total_count:
            print_log(" - could not retrieve lenght for " + str(seqid_total_count - seqinfo.size()) + " sequences")
            seqid_total_count = seqinfo.size()
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        tx = time.time()
        print_log("Extracting taxonomic information from accession2taxid files")
        dowloaded_acc2txid_files = []
        for acc2txid in seq_info_mode:
            if acc2txid not in acc2txid_options:
                print_log(" - " + acc2txid +  " is not a valid option")
            else:
                dowloaded_acc2txid_files.append(get_accession2taxid(acc2txid, tmp_output_folder))
            
        count_acc2txid = seqinfo.parse_acc2txid(dowloaded_acc2txid_files)
        for acc2txid_file, cnt in count_acc2txid.items():
            print_log(" - " + str(cnt) + " entries found in the " + acc2txid_file.split("/")[-1] + " file")
        # Check if retrieved taxids are the same as number of inputs, reset counter
        if seqinfo.size() < seqid_total_count:
            print_log(" - could not retrieve taxid for " + str(seqid_total_count - seqinfo.size()) + " accessions")
            seqid_total_count = seqinfo.size()
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        if use_assembly:
            tx = time.time()
            print_log("Retrieving assembly information from NCBI E-utils")
            seqid_file = tmp_output_folder + "seqids.txt"
            seqinfo.write_seqid_file(seqid_file)
            seqinfo.parse_ncbi_eutils(seqid_file, path_exec['get_len_taxid'], skip_len_taxid=True, get_assembly=True)
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

def set_taxdump_files(cfg, tmp_output_folder):
    if not cfg.taxdump_file:
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(get_taxdump(tmp_output_folder), tmp_output_folder)
    elif cfg.taxdump_file[0].endswith(".tar.gz"):
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(cfg.taxdump_file[0], tmp_output_folder)
    else:
        ncbi_nodes_file = cfg.taxdump_file[0]
        ncbi_names_file = cfg.taxdump_file[1]
        ncbi_merged_file =  cfg.taxdump_file[2] if len(cfg.taxdump_file)==3 else ""

    return ncbi_nodes_file, ncbi_merged_file, ncbi_names_file

def get_taxdump(tmp_output_folder):
    tx = time.time()
    print_log("Downloading taxdump")
    taxdump_file = tmp_output_folder+'taxdump.tar.gz'
    run_wget_taxdump_cmd = 'wget -qO {0} "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"'.format(taxdump_file)
    stdout, stderr = run(run_wget_taxdump_cmd, print_stderr=True)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")
    return taxdump_file

def unpack_taxdump(taxdump_file, tmp_output_folder):
    tx = time.time()
    print_log("Unpacking taxdump")
    unpack_taxdump_cmd = 'tar xf {0} -C "{1}" nodes.dmp merged.dmp names.dmp'.format(taxdump_file, tmp_output_folder)
    stdout, stderr = run(unpack_taxdump_cmd, print_stderr=True)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")
    return tmp_output_folder+'nodes.dmp', tmp_output_folder+'names.dmp', tmp_output_folder+'merged.dmp'

def get_accession2taxid(acc2txid, tmp_output_folder):
    tx = time.time()
    acc2txid_file = acc2txid + ".accession2taxid.gz"
    print_log("Downloading " + acc2txid_file)
    acc2txid_file = tmp_output_folder + acc2txid_file
    run_wget_acc2txid_file_cmd = 'wget -qO {0} "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/{1}.accession2taxid.gz"'.format(acc2txid_file, acc2txid)
    stdout, stderr = run(run_wget_acc2txid_file_cmd, print_stderr=True)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")
    return acc2txid_file

def ibf_size_mb(bin_len, n_bins, max_fp, hash_functions, kmer_size):
    return (math.ceil(-(1/((1-max_fp**(1/float(hash_functions)))**(1/float(hash_functions*(bin_len-kmer_size+1)))-1)))*optimal_bins(n_bins))/8388608

def approx_n_bins(bin_len, overlap_len, groups_len): 
    frag_len=bin_len-overlap_len
    n_bins = sum([math.ceil(math.ceil(l/(frag_len-overlap_len))/(bin_len/(frag_len+overlap_len))) for l in groups_len.values()])
    return n_bins

def optimal_bins(n): 
    #return optimal number of bins for the IBF (multiples of 64)
    return (math.floor(n/64)+1)*64 

def estimate_bin_len_size(cfg, seqinfo, tax, use_assembly):
    # Simulate bins that will be created by taxsbp with many bin lenghts
    # Select the best trade-off on size and n. of bins

    # Generate dict with groups:total_length
    groups_len = {}
    if use_assembly:
        groups_len = seqinfo.seqinfo.groupby('assembly').sum().to_dict()['length']
    elif cfg.rank=="taxid":
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

    print_log(" - Approx. min. size possible: " + str("{0:.2f}".format(min_filter_size)) + "MB")
    # Define max size
    # if none defined or too small, Define the max as 1.5 time size of the min
    if cfg.max_bloom_size is None:
        max_filter_size = min_filter_size*1.5
    elif cfg.max_bloom_size<min_filter_size:
        max_filter_size = min_filter_size*1.5
        print_log(" - --max-bloom-size " + str(cfg.max_bloom_size) + "MB is too small, using max. default (1.5x min.): " + str("{0:.2f}".format(max_filter_size)) + "MB")
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

def validate_input_files(cfg):
    input_files_from_directory = []
    input_files = []

    # get files from directory
    if cfg.input_directory and cfg.input_extension:
        if not os.path.isdir(cfg.input_directory):
            print_log(cfg.input_directory + " is not a valid directory")
        else:
            for file in os.listdir(cfg.input_directory):
                if file.endswith(cfg.input_extension):
                    input_files_from_directory.append(os.path.join(cfg.input_directory, file))
            print_log(str(len(input_files_from_directory)) + " file(s) [" + cfg.input_extension + "] found in " + cfg.input_directory)

    # remove non existent files from input list
    if cfg.input_files: 
        input_files = check_files(cfg.input_files)

    return input_files, input_files_from_directory
