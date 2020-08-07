#!/usr/bin/env python3

import argparse, os, sys, subprocess, time, shlex, shutil, gzip, pickle, math
import pandas as pd
import numpy as np
from io import StringIO
from collections import defaultdict, OrderedDict
import taxsbp.taxsbp

def main(arguments=None):

    version = '0.3.0'
    
    ####################################################################################################
	
    # Arguments from testing
    if arguments is not None: sys.argv=arguments
	
    build_parser = argparse.ArgumentParser(description='Build options', add_help=False)
    
    # Required
    build_group_required = build_parser.add_argument_group('required arguments')
    build_group_required.add_argument('-d', '--db-prefix',      required=True, type=str,                    metavar='db_prefix',        help='Database output prefix (.ibf, .map, .tax, .gnn will be created)')
    build_group_required.add_argument('-i', '--input-files',    required=False, type=str, nargs="*",         metavar='',  help='Input reference sequence fasta files [.gz]')
    
    # Defaults
    build_group_optional = build_parser.add_argument_group('optional arguments')
    build_group_optional.add_argument('-r', '--rank',            type=str,   default='species', metavar='', help='Target taxonomic rank for classification [assembly,taxid,species,genus,...]. Default: species')
    build_group_optional.add_argument('-k', '--kmer-size',       type=int,   default=19,      metavar='', help='The k-mer size for the interleaved bloom filter. Default: 19')
    build_group_optional.add_argument('-n', '--hash-functions',  type=int,   default=3,       metavar='', help='The number of hash functions for the interleaved bloom filter. Default: 3')
    build_group_optional.add_argument('-f', '--max-fp',          type=float, default=0.05,    metavar='', help='Max. false positive rate for k-mer classification. Default: 0.05')
    build_group_optional.add_argument('-m', '--max-bloom-size',  type=int,                    metavar='', help='Approx. maximum filter size in Megabytes (MB). Will estimate best --bin-length based on --kmer-size, --hash-functions and --max-fp  [Mutually exclusive --fixed-bloom-size]')
    build_group_optional.add_argument('-l', '--bin-length',      type=int,                    metavar='', help='Maximum length (in bp) for each bin. Default: auto')
    build_group_optional.add_argument('-t', '--threads',         type=int,   default=2,       metavar='', help='Number of subprocesses/threads to use. Default: 2')
    build_group_optional.add_argument('--fixed-bloom-size',      type=int,                    metavar='', help='Fixed size for filter in Megabytes (MB), will ignore --max-fp [Mutually exclusive --max-bloom-size] ')
    build_group_optional.add_argument('--fragment-length',       type=int,   default=-1,      metavar='', help='Fragment length (in bp). Set to 0 to not fragment sequences. Default: --bin-length - --overlap-length')
    build_group_optional.add_argument('--overlap-length',        type=int,   default=300,     metavar='', help='Fragment overlap length (in bp). Should be bigger than the read length used for classification. Default: 300')
    build_group_optional.add_argument('--seq-info-mode',         type=str, nargs="*", default=["auto"],  metavar='', help='Mode to obtain sequence information. For each sequence entry provided, ganon requires taxonomic and seq. length information. If a small number of sequences is provided (<50000) or when --rank assembly, ganon will automatically obtain data with NCBI E-utils websevices (eutils). Offline mode will download batch files from NCBI Taxonomy and look for taxonomic ids in the order provided. Options: [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot], eutils (force webservices) or auto (uses eutils or [nucl_gb nucl_wgs]). Default: auto [Mutually exclusive --seq-info-file]')
    build_group_optional.add_argument('--seq-info-file',         type=str,                               metavar='', help='Pre-generated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly id]) [Mutually exclusive --seq-info]')
    build_group_optional.add_argument('--taxdump-file',          type=str, nargs="*",                    metavar='', help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
    build_group_optional.add_argument('--input-directory',       type=str,                    metavar='', help='Directory containing input files')
    build_group_optional.add_argument('--input-extension',       type=str,                    metavar='', help='Extension of files to use with --input-directory (provide it without * expansion, e.g. ".fna.gz")')

    # Extra
    build_group_optional.add_argument('--verbose', default=False, action='store_true', help='Verbose mode for ganon-build')
    build_group_optional.add_argument('--ganon-path', type=str, default="", help=argparse.SUPPRESS)
    build_group_optional.add_argument('--n-refs', type=int, help=argparse.SUPPRESS)
    build_group_optional.add_argument('--n-batches', type=int, help=argparse.SUPPRESS)

    ####################################################################################################

    update_parser = argparse.ArgumentParser(description='Update options', add_help=False)

    # Required
    update_group_required = update_parser.add_argument_group('required arguments')
    update_group_required.add_argument('-d', '--db-prefix',         required=True,  type=str,               metavar='db_prefix',        help='Database input prefix (.ibf, .map, .tax, .gnn)')
    update_group_required.add_argument('-i', '--input-files',       required=False, type=str, nargs="*",    metavar='',  help='Input reference sequence fasta files [.gz] to be included to the database. Complete set of updated sequences should be provided when using --update-complete')
    
    # Defaults
    update_group_optional = update_parser.add_argument_group('optional arguments')
    update_group_optional.add_argument('-o', '--output-db-prefix',                  type=str,                               metavar='', help='Output database prefix (.ibf, .map, .tax, .gnn). Default: overwrite current --db-prefix')
    update_group_optional.add_argument('-c', '--update-complete',                             default=False, action='store_true', help='Update adding and removing sequences. Input files should represent the complete updated set of references, not only new sequences.')
    update_group_optional.add_argument('-t', '--threads',                           type=int, default=2,                    metavar='', help='Number of subprocesses/threads to use. Default: 2')
    update_group_optional.add_argument('--seq-info-mode',         type=str, nargs="*", default=["auto"],  metavar='', help='Mode to obtain sequence information. For each sequence entry provided, ganon requires taxonomic and seq. length information. If a small number of sequences is provided (<50000) or when --rank assembly, ganon will automatically obtained data with NCBI E-utils websevices (eutils). Offline mode will download batch files from NCBI Taxonomy and look for taxonomic ids in the order provided. Options: [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot], eutils (force webservices) or auto (uses eutils or [nucl_gb nucl_wgs]). Default: auto [Mutually exclusive --seq-info-file]')
    update_group_optional.add_argument('--seq-info-file',         type=str,                               metavar='', help='Pre-generated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly id]) [Mutually exclusive --seq-info]')
    update_group_optional.add_argument('--taxdump-file',          type=str, nargs="*",                    metavar='', help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
    update_group_optional.add_argument('--input-directory',       type=str,                    metavar='', help='Directory containing input files')
    update_group_optional.add_argument('--input-extension',       type=str,                    metavar='', help='Extension of files to use with --input-directory (provide it without * expansion, e.g. ".fna.gz")')

    # Extra
    update_group_optional.add_argument('--verbose', default=False, action='store_true', help='Verbose mode for ganon-build')
    update_group_optional.add_argument('--ganon-path', type=str, default="", help=argparse.SUPPRESS)
    update_group_optional.add_argument('--n-refs', type=int, help=argparse.SUPPRESS)
    update_group_optional.add_argument('--n-batches', type=int, help=argparse.SUPPRESS)

    ####################################################################################################

    classify_parser = argparse.ArgumentParser(description='Classification options', add_help=False)

    # Required
    classify_group_required = classify_parser.add_argument_group('required arguments')
    classify_group_required.add_argument('-d', '--db-prefix', required=True, nargs="*", type=str, metavar='db_prefix', help='Database input prefix[es]')
    classify_group_required.add_argument('-r', '--single-reads', nargs="*", type=str, metavar='reads.fq[.gz]', help='Multi-fastq[.gz] file[s] to classify')
    classify_group_required.add_argument('-p', '--paired-reads', nargs="*", type=str,  metavar='reads.1.fq[.gz] reads.2.fq[.gz]', help='Multi-fastq[.gz] pairs of file[s] to classify')

    # Defaults
    classify_group_optional = classify_parser.add_argument_group('optional arguments')
    classify_group_optional.add_argument('-c', '--hierarchy-labels', type=str,    nargs="*", help='Hierarchy definition, one for each database input. Can also be a string, but input will be sorted to define order (e.g. 1 1 2 3). Default: 1')
    classify_group_optional.add_argument('-k', '--min-kmers',        type=float,  nargs="*", help='Min. percentage of k-mers matching to consider a read assigned. Single value or one per database (e.g. 0.5 0.7 1 0.25). Default: 0.25 [Mutually exclusive --max-error]')
    classify_group_optional.add_argument('-e', '--max-error',        type=int,    nargs="*", help='Max. number of errors allowed. Single value or one per database (e.g. 3 3 4 0) [Mutually exclusive --min-kmers]')
    classify_group_optional.add_argument('-u', '--max-error-unique', type=int,    nargs="*", help='Max. number of errors allowed for unique assignments after filtering. Matches below this error rate will not be discarded, but assigned to a parent taxonomic level. Single value or one per hierarchy (e.g. 0 1 2). -1 to disable. Default: -1')    
    classify_group_optional.add_argument('-l', '--strata-filter',    type=int,    nargs="*", help='Additional errors allowed (relative to the best match) to filter and select matches. Single value or one per hierarchy (e.g. 0 1 2). -1 to disable filtering. Default: 0')    
    classify_group_optional.add_argument('-f', '--offset',           type=int,               help='Number of k-mers to skip during classification. Can speed up analysis but may reduce recall. (e.g. 1 = all k-mers, 3 = every 3rd k-mer). Default: 2')    
    classify_group_optional.add_argument('-o', '--output-prefix',    type=str,               help='Output prefix for .lca and .rep. Empty to output to STDOUT (only .lca will be printed)')
    classify_group_optional.add_argument('-a', '--output-all',          default=False, action='store_true', help='Output an additional file with all matches (.all). File can be very large.')
    classify_group_optional.add_argument('-n', '--output-unclassified', default=False, action='store_true', help='Output an additional file with unclassified read headers (.unc)')
    classify_group_optional.add_argument('-s', '--output-single',       default=False, action='store_true', help='When using multiple hierarchical levels, output everything in one file instead of one per hierarchy')
    classify_group_optional.add_argument('--ranks', type=str, default=[], nargs="*", help='Ranks to show in the report (.tre). "all" for all identified ranks. empty for default ranks: superkingdom phylum class order family genus species species+ assembly. This file can be re-generated with the ganon report command.')

    classify_group_optional.add_argument('-t', '--threads', type=int, help='Number of subprocesses/threads to use. Default: 3')
    classify_group_optional.add_argument('--n-reads', type=int, help=argparse.SUPPRESS)
    classify_group_optional.add_argument('--n-batches', type=int, help=argparse.SUPPRESS)
    classify_group_optional.add_argument('--verbose', default=False, action='store_true',  help='Verbose mode for ganon-classify')
    classify_group_optional.add_argument('--ganon-path', type=str, default="", help=argparse.SUPPRESS) 

    ####################################################################################################

    report_parser = argparse.ArgumentParser(description='Report options', add_help=False)

    # Required
    report_group_required = report_parser.add_argument_group('required arguments')
    report_group_required.add_argument('-i', '--rep-file',  required=True, type=str, help='{prefix}.rep file output from ganon classify')
    report_group_required.add_argument('-d', '--db-prefix', required=True, type=str, nargs="*", metavar='db_prefix', help='Database prefix[es] used for classification.')
    
    # Defaults
    report_group_optional = report_parser.add_argument_group('optional arguments')
    report_group_optional.add_argument('-r', '--ranks', type=str, default=[], nargs="*", help='Ranks for the final report. "all" for all identified ranks. empty for default ranks: superkingdom phylum class order family genus species species+ assembly')
    report_group_optional.add_argument('-m', '--min-matches', type=int, default=0, help='Min. number of matches to output. 0 for all. Default: 0')
    report_group_optional.add_argument('-p', '--min-matches-perc', type=float, default=0, help='Min. percentage of matches to output. 0 for all. Default: 0')
    report_group_optional.add_argument('-t', '--taxids', type=str, default=[], nargs="*", help='One or more taxids to filter report. Example: 562 2157 report only E. Coli and Archaea matches')
    report_group_optional.add_argument('-o', '--output-report', type=str, help='Output file for report. Default: STDOUT')

    ####################################################################################################

    parser = argparse.ArgumentParser(prog='ganon', description='ganon', conflict_handler='resolve')
    parser.add_argument('-v', '--version', action='version', version='version: %(prog)s ' + version, help="Show program's version number and exit.")
    subparsers = parser.add_subparsers()
    
    build = subparsers.add_parser('build', help='Build ganon database', parents=[build_parser])
    build.set_defaults(which='build')

    update = subparsers.add_parser('update', help='Update ganon database', parents=[update_parser])
    update.set_defaults(which='update')

    classify = subparsers.add_parser('classify', help='Classify reads', parents=[classify_parser])
    classify.set_defaults(which='classify')

    report = subparsers.add_parser('report', help='Generate reports', parents=[report_parser])
    report.set_defaults(which='report')

    args = parser.parse_args()

    if len(sys.argv[1:])==0: # Print help calling script without parameters
        parser.print_help() 
        sys.exit(0)
    
    # set path for executables
    path_exec = set_paths(args)
    if not path_exec: sys.exit(1)
    # validate arguments
    if not validate_args(args): sys.exit(1)

    if args.which in ['build','update']:    
        # validate input files
        input_files, input_files_from_directory = validate_input_files(args)
        if len(input_files)==0 and len(input_files_from_directory)==0:
            print_log("No valid input files found")
            sys.exit(1)

        # set output files
        db_prefix = args.db_prefix
        if args.which=='update' and args.output_db_prefix:
            tmp_output_folder = args.output_db_prefix + "_tmp/"
        else:
            tmp_output_folder = db_prefix + "_tmp/"
        db_prefix_ibf = db_prefix + ".ibf"
        db_prefix_map = db_prefix + ".map"
        db_prefix_tax = db_prefix + ".tax"
        db_prefix_gnn = db_prefix + ".gnn"


    tx_total = time.time()
    print_log("- - - - - - - - - -")
    print_log("   _  _  _  _  _   ")
    print_log("  (_|(_|| |(_)| |  ")
    print_log("   _|   v. "+ str(version))
    print_log("- - - - - - - - - -")

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    if args.which=='build':              
        set_tmp_folder(tmp_output_folder)

        # Set assembly mode
        use_assembly=True if args.rank=="assembly" else False

        # Set up taxonomy
        ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(args, tmp_output_folder)
        tx = time.time()
        print_log("Parsing taxonomy")
        tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        # Load seqids and generate seqinfo
        if args.seq_info_file:
            seqinfo = load_seqids(seq_info_file=args.seq_info_file)
        else:
            seqinfo = load_seqids(files=input_files + input_files_from_directory) 
            load_seqinfo(tmp_output_folder, seqinfo, path_exec, args.seq_info_mode, use_assembly)

        # Set bin length
        if args.bin_length: # user defined
            bin_length = args.bin_length
        else:
            tx = time.time()
            print_log("Calculating best bin length")
            bin_length = estimate_bin_len2(args, seqinfo, tax, use_assembly)
            if bin_length<=0: 
                bin_length=1000000
                print_log(" - could not estimate bin length, using default (" + str(bin_length) + ")")
            else:
                print_log(" - " + str(bin_length) + "bp")
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        # Set fragment length
        if args.fragment_length==-1: # if ==-1 set default
            fragment_length = bin_length - args.overlap_length
        elif args.fragment_length==0: # if ==0 deactivate
            fragment_length = 0
        else: # user input
            fragment_length = args.fragment_length - args.overlap_length 

        tx = time.time()
        print_log("Running taxonomic clustering (TaxSBP)")
        taxsbp_params={}
        taxsbp_params["nodes_file"] = ncbi_nodes_file
        taxsbp_params["bin_len"] = bin_length
        if use_assembly:
            taxsbp_params["bin_exclusive"] = "assembly"
        elif args.rank=="taxid":
            taxsbp_params["bin_exclusive"] = "leaves"
        else:
            taxsbp_params["bin_exclusive"] =  args.rank
        if ncbi_merged_file: taxsbp_params["merged_file"] = ncbi_merged_file
        if fragment_length: 
            taxsbp_params["fragment_len"] = fragment_length
            taxsbp_params["overlap_len"] = args.overlap_length
        if use_assembly: taxsbp_params["specialization"] = "assembly"
        taxsbp_params["input_table"] = seqinfo.get_csv()
        bins = Bins(taxsbp_ret=taxsbp.taxsbp.pack(**taxsbp_params))
        del taxsbp_params
        # bin statistics
        actual_number_of_bins = bins.get_number_of_bins()
        optimal_number_of_bins = optimal_bins(actual_number_of_bins)
        max_length_bin = bins.get_max_bin_length()
        max_kmer_count = max_length_bin-args.kmer_size+1 # aproximate number of unique k-mers by just considering that they are all unique
        print_log(" - " + str(actual_number_of_bins) + " bins created")
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        tx = time.time()
        print_log("Building database files")
        
        # Write .map file
        print_log(" - " + db_prefix_map)
        bins.write_map_file(db_prefix_map, use_assembly)

        # Write .tax file
        print_log(" - " + db_prefix_tax)
        tax.filter(bins.get_taxids()) # filter only used taxids
        if use_assembly: tax.add_nodes(bins.get_specialization_taxid(), "assembly") # add assembly nodes
        tax.write(db_prefix_tax)

        # Write .gnn file
        print_log(" - " + db_prefix_gnn)
        gnn = Gnn(kmer_size=args.kmer_size, 
                hash_functions=args.hash_functions, 
                number_of_bins=actual_number_of_bins, 
                rank=args.rank,
                bin_length=bin_length,
                fragment_length=fragment_length,
                overlap_length=args.overlap_length,
                bins=bins.get_list())
        gnn.write(db_prefix_gnn)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        print_log("Building index (ganon-build)")
        # define bloom filter size based on given false positive
        MBinBits = 8388608
        print_log(" - max unique " + str(args.kmer_size) +  "-mers: " + str(max_kmer_count))
        if not args.fixed_bloom_size:
            bin_size_bits = math.ceil(-(1/((1-args.max_fp**(1/float(args.hash_functions)))**(1/float(args.hash_functions*max_kmer_count))-1)))   
            print_log(" - IBF calculated size with fp<=" + str(args.max_fp) + ": " + str("{0:.4f}".format((bin_size_bits*optimal_number_of_bins)/MBinBits)) + "MB (" + str(bin_size_bits) + " bits/bin * " + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])")
        else:
            bin_size_bits = math.ceil((args.fixed_bloom_size * MBinBits)/optimal_number_of_bins);
            estimated_max_fp = (1-((1-(1/float(bin_size_bits)))**(args.hash_functions*max_kmer_count)))**args.hash_functions
            print_log(" - IBF calculated max. fp with size=" + str(args.fixed_bloom_size) + "MB: " + str("{0:.4f}".format(estimated_max_fp) + " ("  + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])"))

        # Write aux. file for ganon
        acc_bin_file = tmp_output_folder + "acc_bin.txt"
        bins.write_acc_bin_file(acc_bin_file)

        run_ganon_build_cmd = " ".join([path_exec['build'],
                                        "--seqid-bin-file " + acc_bin_file,
                                        "--filter-size-bits " + str(bin_size_bits*optimal_number_of_bins) if args.max_fp else "--filter-size " + str(args.fixed_bloom_size),
                                        "--kmer-size " + str(args.kmer_size),
                                        "--hash-functions " + str(args.hash_functions),
                                        "--threads " + str(args.threads),
                                        "--output-filter-file " + db_prefix_ibf,
                                        "--verbose" if args.verbose else "",
                                        "--n-refs " + str(args.n_refs) if args.n_refs is not None else "",
                                        "--n-batches " + str(args.n_batches) if args.n_batches is not None else "",
                                        "--reference-files " + ",".join([file for file in input_files]) if input_files else "",
                                        "--directory-reference-files " + args.input_directory if args.input_directory else "",
                                        "--extension " + args.input_extension if args.input_extension else ""])
        stdout, stderr = run(run_ganon_build_cmd, print_stderr=True)

        # Delete temp files
        shutil.rmtree(tmp_output_folder)
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    elif args.which=='update':  
        tx = time.time()
        set_tmp_folder(tmp_output_folder)

        # Load .gnn file   
        gnn = Gnn(file=db_prefix_gnn)
        # Set assembly mode
        use_assembly=True if gnn.rank=="assembly" else False

        # load bins
        bins = Bins(taxsbp_ret=gnn.bins)

        # Load seqids and generate seqinfo
        if args.seq_info_file:
            seqinfo = load_seqids(seq_info_file=args.seq_info_file)
        else:
            seqinfo = load_seqids(files=input_files + input_files_from_directory) 

        # check sequences compared to bins
        added_seqids, removed_seqids, kept_seqids = check_updated_seqids(set(seqinfo.get_seqids()), set(bins.get_seqids()))
        # Ignore removed sequences if not doing complete update
        if args.update_complete: 
            print_log("Update: adding " + str(len(added_seqids)) + " sequences, removing " + str(len(removed_seqids)) + " sequences, keeping " + str(len(kept_seqids)) + " sequences")
        else:
            removed_seqids=[]
            print_log("Update: adding " + str(len(added_seqids)) + " sequences, ignoring " + str(len(kept_seqids)) + " repeated sequences")
        print_log("")

        if not added_seqids and not removed_seqids:
            print_log("Nothing to update.")
            sys.exit(0)

        if args.update_complete:           
            # Remove already included seqids to just retrieve information for added sequences
            seqinfo.remove_seqids(kept_seqids | removed_seqids)
        else:
            # Remove seqids already present in the current version (repeated entries)
            seqinfo.remove_seqids(kept_seqids)

        # load seqinfo file with data (after removing ids)
        if not args.seq_info_file: load_seqinfo(tmp_output_folder, seqinfo, path_exec, args.seq_info_mode, use_assembly)

        # save set of current binids
        previous_binids = set(bins.get_binids())
        # remove seqids from bins if performing update complete
        if args.update_complete and removed_seqids:
            bins.remove_seqids(removed_seqids)
        # save set of kept binids after removal
        kept_binids = set(bins.get_binids())

        # Set up taxonomy files
        ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = set_taxdump_files(args, tmp_output_folder)

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
        print_log(" - " + args.output_db_prefix + ".tax" if args.output_db_prefix else db_prefix_tax)
        tax = Tax(ncbi_nodes=ncbi_nodes_file, ncbi_names=ncbi_names_file)
        # Update and write .tax file
        tax.filter(updated_bins.get_taxids()) # filter only used taxids
        if use_assembly: tax.add_nodes(updated_bins.get_specialization_taxid(), "assembly") # add assembly nodes
        # Load old .tax file into new taxonomy
        tax.merge(Tax([db_prefix_tax]))
        # Write .tax file
        tax.write(args.output_db_prefix + ".tax" if args.output_db_prefix else db_prefix_tax)
        # TODO - remove entries from .tax from removed entries of the db

        # merge updated and old bins together
        bins.merge(updated_bins)
        
        # Write .gnn file
        print_log(" - " + args.output_db_prefix + ".gnn" if args.output_db_prefix else db_prefix_gnn)
        gnn.bins = bins.get_list() # save updated bins
        gnn.number_of_bins=bins.get_number_of_bins() # add new bins count
        gnn.write(args.output_db_prefix + ".gnn" if args.output_db_prefix else db_prefix_gnn)

        # Recreate .map file based on the new bins
        print_log(" - " + args.output_db_prefix + ".map" if args.output_db_prefix else db_prefix_map)
        bins.write_map_file(args.output_db_prefix + ".map" if args.output_db_prefix else db_prefix_map, use_assembly)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")

        tx = time.time()
        print_log("Updating index (ganon-build)")

        # Write aux. file for ganon
        # This file has to contain all new sequences
        # in case of update_complete, 
        acc_bin_file = tmp_output_folder + "acc_bin.txt"
        
        if args.update_complete:
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
        run_ganon_build_cmd = " ".join([path_exec['build'],
                                        "--update-filter-file " + db_prefix_ibf,
                                        "--seqid-bin-file " + acc_bin_file,
                                        "--output-filter-file " + tmp_db_prefix_ibf,
                                        "--threads " + str(args.threads),                                
                                        "--verbose" if args.verbose else "",
                                        "--n-refs " + str(args.n_refs) if args.n_refs is not None else "",
                                        "--n-batches " + str(args.n_batches) if args.n_batches is not None else "",
                                        "--reference-files " + ",".join([file for file in input_files]) if input_files else "",
                                        "--directory-reference-files " + args.input_directory if args.input_directory else "",
                                        "--extension " + args.input_extension if args.input_extension else "",
                                        "--update-complete" if args.update_complete else ""])
        stdout, stderr = run(run_ganon_build_cmd, print_stderr=True)
        
        # move IBF to final location
        shutil.move(tmp_db_prefix_ibf, args.output_db_prefix + ".ibf" if args.output_db_prefix else db_prefix_ibf)

        # Delete temp files
        shutil.rmtree(tmp_output_folder)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    elif args.which=='classify':

        print_log("Classifying reads (ganon-classify)")
        run_ganon_classify = " ".join([path_exec['classify'],
                                       "--single-reads " +  ",".join(args.single_reads) if args.single_reads else "",
                                       "--paired-reads " +  ",".join(args.paired_reads) if args.paired_reads else "",
                                       "--ibf " + ",".join([db_prefix+".ibf" for db_prefix in args.db_prefix]),
                                       "--map " + ",".join([db_prefix+".map" for db_prefix in args.db_prefix]), 
                                       "--tax " + ",".join([db_prefix+".tax" for db_prefix in args.db_prefix]),
                                       "--hierarchy-labels " + ",".join(args.hierarchy_labels) if args.hierarchy_labels else "",
                                       "--max-error " + ",".join([str(me) for me in args.max_error]) if args.max_error else "",
                                       "--min-kmers " + ",".join([str(mk) for mk in args.min_kmers]) if args.min_kmers else "",
                                       "--max-error-unique " + ",".join([str(meu) for meu in args.max_error_unique]) if args.max_error_unique else "",
                                       "--strata-filter " + ",".join([str(sf) for sf in args.strata_filter]) if args.strata_filter else "",
                                       "--offset " + str(args.offset) if args.offset else "",
                                       "--output-prefix " + args.output_prefix if args.output_prefix else "",
                                       "--output-all" if args.output_all else "",
                                       "--output-unclassified" if args.output_unclassified else "",
                                       "--output-single" if args.output_single else "",
                                       "--threads " + str(args.threads) if args.threads else "",
                                       "--n-reads " + str(args.n_reads) if args.n_reads is not None else "",
                                       "--n-batches " + str(args.n_batches) if args.n_batches is not None else "",
                                       "--verbose" if args.verbose else ""])
        stdout, stderr = run(run_ganon_classify)
        if not args.output_prefix: print(stdout)
        print_log(stderr)

        if args.output_prefix:
            tx = time.time()
            print_log("Generating report")
            tax = Tax([db_prefix+".tax" for db_prefix in args.db_prefix])
            classified_reads, unclassified_reads, reports = parse_rep(args.output_prefix+".rep")
            print_final_report(reports, tax, classified_reads, unclassified_reads, args.output_prefix+".tre", args.ranks, 0, 0, [])
            print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n")


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    elif args.which=='report':
        classified_reads, unclassified_reads, reports = parse_rep(args.rep_file)
        tax = Tax([db_prefix+".tax" for db_prefix in args.db_prefix])
        print_final_report(reports, tax, classified_reads, unclassified_reads, args.output_report, args.ranks, args.min_matches, args.min_matches_perc, args.taxids)

    print_log("Total elapsed time: " + str("%.2f" % (time.time() - tx_total)) + " seconds.")

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def run(cmd, output_file=None, print_stderr=False, shell=False, exit_on_error: bool=True):
    errcode=0
    stdout=""
    stderr=""
    try:
        process = subprocess.Popen(shlex.split(cmd) if not shell else cmd, 
                                    shell=shell, 
                                    universal_newlines=True, 
                                    stdout=subprocess.PIPE if output_file is None else open(output_file, 'w'), 
                                    stderr=subprocess.PIPE)   
        stdout, stderr = process.communicate() # wait for the process to terminate
        errcode = process.returncode
        if errcode!=0: raise Exception()
        if print_stderr and stderr: print_log(stderr)
 
    #except OSError as e: # The most common exception raised is OSError. This occurs, for example, when trying to execute a non-existent file. Applications should prepare for OSError exceptions.
    #except ValueError as e: #A ValueError will be raised if Popen is called with invalid arguments.
    except Exception as e:
        print_log('The following command failed to execute:\n'+cmd)
        print_log(str(e))
        print_log("Error code: "+str(errcode))
        print_log("Out: ")
        if stdout: print_log(stdout)
        print_log("Error: ")
        if stderr: print_log(stderr)
        if exit_on_error: sys.exit(errcode)

    return stdout, stderr

def check_updated_seqids(new_seqids, old_seqids):
    # remove repeated from old bins
    added_seqids = new_seqids.difference(old_seqids)
    removed_seqids = old_seqids.difference(new_seqids)
    kept_seqids = old_seqids.difference(removed_seqids)

    return added_seqids, removed_seqids, kept_seqids

def set_tmp_folder(tmp_output_folder):
    # Create temporary working directory
    if os.path.exists(tmp_output_folder): shutil.rmtree(tmp_output_folder) # delete if already exists
    os.makedirs(tmp_output_folder)

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

def set_taxdump_files(args, tmp_output_folder):
    if not args.taxdump_file:
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(get_taxdump(tmp_output_folder), tmp_output_folder)
    elif args.taxdump_file[0].endswith(".tar.gz"):
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(args.taxdump_file[0], tmp_output_folder)
    else:
        ncbi_nodes_file = args.taxdump_file[0]
        ncbi_names_file = args.taxdump_file[1]
        ncbi_merged_file =  args.taxdump_file[2] if len(args.taxdump_file)==3 else ""

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

def parse_rep(rep_file):
    reports = defaultdict(lambda: defaultdict(lambda: {'direct_matches': 0, 'unique_reads': 0, 'lca_reads': 0}))
    with open(rep_file, 'r') as rep_file:
        for line in rep_file:
            fields = line.rstrip().split("\t")
            if fields[0] == "#total_classified":
                seq_cla = int(fields[1])
            elif fields[0] == "#total_unclassified":
                seq_unc = int(fields[1])
            else:
                hierarchy_name, target, direct_matches, unique_reads, lca_reads, rank, name = fields
                reports[hierarchy_name][target]["direct_matches"]+=int(direct_matches)
                reports[hierarchy_name][target]["unique_reads"]+=int(unique_reads)
                reports[hierarchy_name][target]["lca_reads"]+=int(lca_reads)
    return seq_cla, seq_unc, reports

def ibf_size_mb(bin_len, n_bins, max_fp, hash_functions, kmer_size):
    return (math.ceil(-(1/((1-max_fp**(1/float(hash_functions)))**(1/float(hash_functions*(bin_len-kmer_size+1)))-1)))*n_bins)/8388608

def approx_n_bins(bin_len, overlap_len, groups_len): 
    frag_len=bin_len-overlap_len
    n_bins = sum([math.ceil(math.ceil(l/(frag_len-overlap_len))/(bin_len/(frag_len+overlap_len))) for l in groups_len.values()])
    return (math.floor(n_bins/64)+1)*64 #return optimal number of bins for the IBF (multiples of 64)

def estimate_bin_len2(args, seqinfo, tax, use_assembly):
    groups_len = {}
    if use_assembly:
        groups_len = seqinfo.seqinfo.groupby('assembly').sum().to_dict()['length']
    else:
        groups_len = pd.concat([seqinfo.seqinfo['taxid'].apply(lambda x: tax.get_rank(x, args.rank)), seqinfo.seqinfo['length']], axis=1).groupby('taxid').sum().to_dict()['length']

    # Set limits
    min_bin_len = min(groups_len.values())+args.overlap_length
    max_bin_len = max(groups_len.values())
    
    # Define intervals of bin_length
    # geometric space to have more simulations on smaller filters
    bin_lens = np.geomspace(min_bin_len, max_bin_len, num=200)
    filter_sizes = np.array([ibf_size_mb(b, approx_n_bins(b, args.overlap_length, groups_len), args.max_fp, args.hash_functions, args.kmer_size) for b in bin_lens])

    # keep only positives
    pos_idx = filter_sizes>0
    bin_lens = bin_lens[pos_idx]
    filter_sizes = filter_sizes[pos_idx]
    if args.max_bloom_size:
        idx_below_max_size = filter_sizes<=args.max_bloom_size
        # if maxsize chosen is smaller than max, recalculate curve
        if idx_below_max_size.any() and sum(idx_below_max_size)>1: 
            # linear space of subset
            bin_lens = np.linspace(bin_lens[idx_below_max_size].min(), bin_lens[idx_below_max_size].max(), 200)
            filter_sizes = np.array([ibf_size_mb(b, approx_n_bins(b, args.overlap_length, groups_len), args.max_fp, args.hash_functions, args.kmer_size) for b in bin_lens])
        else:
            print("max size too small, min around " + str(min(filter_sizes)) + "MB")

    return int(np.percentile(bin_lens, 75))

def print_final_report(reports, tax, classified_reads, unclassified_reads, final_report_file, ranks, min_matches, min_matches_perc, taxids):
    if not ranks:  
        all_ranks = False
        fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+','assembly']
    elif ranks[0]=="all":
        all_ranks = True
        fixed_ranks = []
    else:
        all_ranks = False
        fixed_ranks = ['root'] + ranks

    # sum counts of each report
    merged_rep = defaultdict(int)
    for rep in reports.values():
        for leaf in rep.keys():
            reads_assigned = rep[leaf]['unique_reads']+rep[leaf]['lca_reads']
            if reads_assigned>0:
                merged_rep[leaf]+=reads_assigned

    final_rep = defaultdict(lambda: {'count': 0, 'rank': ""})
    
    # make cummulative sum of the counts on the lineage
    for leaf in merged_rep.keys():
        count = merged_rep[leaf]
        if all_ranks: # use all nodes on the tree
            t = leaf
            r = tax.nodes[t][1]
        else: # use only nodes of the fixed ranks
            t, r = tax.get_node_rank_fixed(leaf, fixed_ranks)

        while t!="0":
            final_rep[t]['count']+=count
            final_rep[t]['rank']=r
            if all_ranks:
                t = tax.nodes[t][0]
                r = tax.nodes[t][1] if t!="0" else ""
            else:
                t, r = tax.get_node_rank_fixed(tax.nodes[t][0], fixed_ranks)

    # build lineage after all entries were defined
    lineage = {}
    for assignment in final_rep.keys():
        lineage[assignment]=[]
        if all_ranks:
            t=assignment
            while t!="0":
                lineage[assignment].insert(0,t)
                t = tax.nodes[t][0]
        else:
            t, r = tax.get_node_rank_fixed(assignment, fixed_ranks)
            max_rank_idx = fixed_ranks.index(r) # get index of current rank
            while t!="0":
                # Add empty || if fixed rank is missing
                for i in range(max_rank_idx-fixed_ranks.index(r)):
                    lineage[assignment].insert(0,"")
                    max_rank_idx-=1

                lineage[assignment].insert(0,t)
                max_rank_idx-=1
                t, r = tax.get_node_rank_fixed(tax.nodes[t][0], fixed_ranks)

        # if taxids is provided, just keep entries with them (and root)
        if taxids and assignment!="1":
            if not any(t in taxids for t in lineage[assignment]):
                del lineage[assignment]

    total_reads = classified_reads + unclassified_reads
    frfile = open(final_report_file, 'w') if final_report_file else None
    print("unclassified" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ str(unclassified_reads) +"\t"+ str("%.5f" % ((unclassified_reads/total_reads)*100)), file=frfile)
    
    
    if all_ranks:
        sorted_assignments = sorted(lineage, key=lineage.get)
    else:
        sorted_assignments = sorted(lineage, key=lambda k: (fixed_ranks.index(final_rep[k]['rank']), -final_rep[k]['count']), reverse=False)
    
    for assignment in sorted_assignments:
        rank=final_rep[assignment]['rank'] 
        name=tax.nodes[assignment][2]
        reads_unique = rep[assignment]['unique_reads']
        reads = merged_rep[assignment] 
        matches=final_rep[assignment]['count']
        if matches < min_matches: continue
        matches_perc=(final_rep[assignment]['count']/total_reads)*100
        if matches_perc < min_matches_perc: continue
        print(rank, assignment, "|".join(lineage[assignment]), name, reads_unique, reads, matches, "%.5f" % matches_perc, file=frfile, sep="\t")
    
    if final_report_file: frfile.close()

def set_paths(args):
    path_exec = {'build': "", 'classify': "", 'get_len_taxid': ""}

    if args.which=='build' or args.which=='update':
        args.ganon_path = args.ganon_path + "/" if args.ganon_path else ""

        # if path is given, look for binaries only there
        ganon_build_paths = [args.ganon_path, args.ganon_path+"build/"] if args.ganon_path else [None, "build/"]
        for p in ganon_build_paths:
            path_exec['build'] = shutil.which("ganon-build", path=p)
            if path_exec['build'] is not None: break
        if path_exec['build'] is None:
            print_log("ganon-build binary was not found. Please inform a specific path with --ganon-path")
            return False

        ganon_get_len_taxid_paths = [args.ganon_path, args.ganon_path+"scripts/", args.ganon_path+"../scripts/"] if args.ganon_path else [None, "scripts/"]
        for p in ganon_get_len_taxid_paths:
            path_exec['get_len_taxid'] = shutil.which("ganon-get-len-taxid.sh", path=p)
            if path_exec['get_len_taxid'] is not None: break
        if path_exec['get_len_taxid'] is None:
            print_log("ganon-get-len-taxid.sh script was not found. Please inform a specific path with --ganon-path")
            return False

    elif args.which=='classify':
        args.ganon_path = args.ganon_path + "/" if args.ganon_path else ""

        ganon_classify_paths = [args.ganon_path, args.ganon_path+"build/"] if args.ganon_path else [None, "build/"]
        for p in ganon_classify_paths:
            path_exec['classify'] = shutil.which("ganon-classify", path=p)
            if path_exec['classify'] is not None: break
        if path_exec['classify'] is None:
            print_log("ganon-classify binary was not found. Please inform a specific path with --ganon-path")
            return False

    return path_exec

def validate_input_files(args):
    input_files_from_directory = []
    input_files = []

    # get files from directory
    if args.input_directory and args.input_extension:
        if not os.path.isdir(args.input_directory):
            print_log(args.input_directory + " is not a valid directory")
        else:
            for file in os.listdir(args.input_directory):
                if file.endswith(args.input_extension):
                    input_files_from_directory.append(os.path.join(args.input_directory, file))
            print_log(str(len(input_files_from_directory)) + " file(s) [" + args.input_extension + "] found in " + args.input_directory)

    # remove non existent files from input list
    if args.input_files: 
        input_files = check_files(args.input_files)

    return input_files, input_files_from_directory

def validate_args(args):
    if args.which in ['build','update']:
        if args.taxdump_file and ((len(args.taxdump_file)==1 and not args.taxdump_file[0].endswith(".tar.gz")) or len(args.taxdump_file)>3):
            print_log("Please provide --taxdump-file taxdump.tar.gz or --taxdump-file nodes.dmp names.dmp [merged.dmp] or leave it empty for automatic download")
            return False

        if not args.input_files and not args.input_directory:
            print_log("Please provide files with --input-files and/or --input-directory with --input-extension")
            return False
        elif args.input_directory and not args.input_extension:
            print_log("Please provide the --input-extension when using --input-directory")
            return False
        elif args.input_directory and "*" in args.input_extension:
            print_log("Please do not use wildcards (*) in the --input-extension")
            return False

        if args.which=='update':
            if not check_db(args.db_prefix):
                return False

        if args.which=='build': #If set (!=0), should be smaller than fragment
            if args.fragment_length>0 and args.overlap_length > args.fragment_length:
                print_log("--overlap-length cannot be bigger than --fragment-length")
                return False

            if args.fixed_bloom_size and not args.bin_length:
                print_log("please set the --bin-length to use --fixed-bloom-size")
                return False

            if args.max_fp<=0:
                print_log("--max-fp has to be bigger than 0")
                return False
        
    elif args.which=='classify':
        for prefix in args.db_prefix:
            if not check_db(prefix):
                return False

        if not args.single_reads and not args.paired_reads:
            print_log("Please provide file[s] with --single-reads or --paired-reads")
            return False

        len_single_reads = 0
        if args.single_reads: 
            args.single_reads = check_files(args.single_reads)
            len_single_reads = len(args.single_reads)
        len_paired_reads = 0
        if args.paired_reads: 
            args.paired_reads = check_files(args.paired_reads)
            len_paired_reads = len(args.paired_reads)
        
        if len_single_reads+len_paired_reads==0:
            print_log("No valid input files found")
            return False

    elif args.which=='report':
        for prefix in args.db_prefix:
            if not check_db(prefix):
                return False

        if not os.path.isfile(args.rep_file):
            print_log("File not found [" + args.rep_file + "]")
            return False

    return True

def check_files(files):
    checked_files = [file for file in files if os.path.isfile(file)]
    if len(checked_files)<len(files):
        print_log(str(len(files)-len(checked_files)) + " input file(s) could not be found")
    return checked_files

def check_db(prefix):
    for db_file_type in [".ibf", ".map", ".tax", ".gnn"]:
        if not os.path.isfile(prefix+db_file_type):
            print_log("Incomplete database [" + prefix  + "] (.ibf, .map, .tax and .gnn)")
            return False
    return True

def print_log(text):
    sys.stderr.write(text+"\n")
    sys.stderr.flush()

class Gnn:
    def __init__(self,
                kmer_size: int = None, 
                hash_functions: int = None, 
                number_of_bins: int = None,
                bin_length: int = None,
                fragment_length: int = None,
                overlap_length: int = None,
                rank: str = None,
                bins: list = None,
                file: str=None):
        if file:
            self.parse(file)
        else:
            self.kmer_size = kmer_size
            self.hash_functions = hash_functions
            self.number_of_bins = number_of_bins
            self.rank = rank
            self.bin_length = bin_length
            self.fragment_length = fragment_length
            self.overlap_length = overlap_length
            self.bins = bins

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'Gnn({})'.format(', '.join(args))

    def parse(self, file):
        gnn = pickle.load(gzip.open(file, "rb"))
        self.kmer_size = gnn['kmer_size']
        self.hash_functions = gnn['hash_functions']
        self.number_of_bins = gnn['number_of_bins']
        self.rank = gnn['rank']
        self.bin_length = gnn['bin_length']
        self.fragment_length = gnn['fragment_length']
        self.overlap_length = gnn['overlap_length']
        self.bins = gnn['bins']

    def write(self, file):
        gnn = {'kmer_size': self.kmer_size, 
        'hash_functions': self.hash_functions, 
        'number_of_bins': self.number_of_bins, 
        'rank': self.rank,
        'bin_length': self.bin_length,
        'fragment_length': self.fragment_length,
        'overlap_length': self.overlap_length,
        'bins': self.bins }
        with gzip.open(file, 'wb') as f: pickle.dump(gnn, f)

class SeqInfo:
    seq_info_colums=['seqid','length','taxid', 'assembly']
    seq_info_types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'assembly': 'str'}
    seqinfo = pd.DataFrame(columns=seq_info_colums)

    def __init__(self, seq_info_file: str=None):
        if seq_info_file:
            self.seqinfo = pd.read_csv(seq_info_file, sep='\t', header=None, skiprows=0, names=self.seq_info_colums, dtype=self.seq_info_types)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'SeqInfo({})'.format(', '.join(args))
    
    def get_csv(self):
        return self.seqinfo.to_csv(sep="\t",header=False, index=False)

    def size(self):
        return self.seqinfo.shape[0]

    def get_seqids(self):
        return self.seqinfo.seqid

    def remove_seqids(self, seqids):
        self.seqinfo = self.seqinfo[~self.seqinfo['seqid'].isin(seqids)]

    def write_seqid_file(self, seqid_file):
        self.seqinfo["seqid"].to_csv(seqid_file, header=False, index=False)

    def parse_seqid(self, input_files):
        for file in input_files:
            # cat | zcat | gawk -> compability with osx
            run_get_header = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{print substr($1,2)}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get_header, print_stderr=False, shell=True)
            self.seqinfo = self.seqinfo.append(pd.read_csv(StringIO(stdout), header=None, names=['seqid']), ignore_index=True)

    def parse_seqid_length(self, input_files):
        for file in input_files:
            # cat | zcat | gawk -> compability with osx
            run_get_length = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{if (seqlen){{print seqlen}}; printf substr($1,2)\"\\t\";seqlen=0;next;}} {{seqlen+=length($0)}}END{{print seqlen}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get_length, print_stderr=False, shell=True)
            self.seqinfo = self.seqinfo.append(pd.read_csv(StringIO(stdout), sep="\t", header=None, names=['seqid', 'length']), ignore_index=True)
            # remove entries with length 0
            self.seqinfo = self.seqinfo[self.seqinfo['length']>0]

    def parse_ncbi_eutils(self, seqid_file, path_exec_get_len_taxid, skip_len_taxid=False, get_assembly=False):
        run_get_len_taxid_cmd = '{0} -i {1} {2} {3} -r'.format(
                                    path_exec_get_len_taxid,
                                    seqid_file,
                                    "-a" if get_assembly else "",
                                    "-s" if skip_len_taxid else "")
        stdout, stderr = run(run_get_len_taxid_cmd, print_stderr=True, exit_on_error=False)
        
        if get_assembly and skip_len_taxid:
            # todo - order is always the same?
            self.seqinfo.assembly = pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=['seqid','assembly'], dtype=self.seq_info_types)['assembly']
        else:
            self.seqinfo = pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=self.seq_info_colums, dtype=self.seq_info_types)

    def parse_acc2txid(self, acc2txid_files):
        count_acc2txid = {}
        set_seqids = set(self.seqinfo.seqid)
        for acc2txid in acc2txid_files:
            tmp_seqid_taxids = pd.read_csv(acc2txid, sep='\t', header=None, skiprows=1, usecols=[1,2], names=['seqid','taxid'], converters={'seqid':lambda x: x if x in set_seqids else ""}, dtype={'taxid': 'str'})
            tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['seqid']!=""] #keep only seqids used
            tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['taxid']!="0"] # filter out taxid==0
            # save count to return
            count_acc2txid[acc2txid] = tmp_seqid_taxids.shape[0]
            # append to main file (will match col names)
            self.seqinfo = self.seqinfo.append(tmp_seqid_taxids, sort=False, ignore_index=True)
            del tmp_seqid_taxids
            #if already found all taxid for the seqids (no need to parse all files till the end)
            if sum(count_acc2txid.values()) == self.size(): 
                break 
        return count_acc2txid


class Bins:
    # bins columns pandas dataframe
    columns=['seqid', 'seqstart', 'seqend', 'length', 'taxid', 'binid', 'specialization']
    bins = pd.DataFrame([], columns=columns)

    def __init__(self, taxsbp_ret: list=[]):
        if taxsbp_ret: 
            self.parse_bins(taxsbp_ret)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'Bins({})'.format(', '.join(args))
    
    def parse_bins(self, taxsbp_ret):
        self.bins = pd.DataFrame(taxsbp_ret, columns=self.columns)

    def size(self):
        return self.bins.shape[0]

    def get_subset(self,seqids: list=[], binids: list=[]):
        subBins = Bins()
        if seqids: subBins.bins = self.bins.loc[self.bins['seqid'].isin(seqids)]
        elif binids: subBins.bins = self.bins.loc[self.bins['binid'].isin(binids)]
        return subBins

    def get_csv(self):
        return self.bins.to_csv(sep="\t",header=False, index=False)

    def get_list(self):
        return self.bins.values.tolist()

    def get_bins(self):
        return self.bins

    def remove_seqids(self, seqids):
        self.bins = self.bins.loc[~self.bins['seqid'].isin(seqids)]

    def get_number_of_bins(self):
        return self.get_binids().size
    
    def get_taxids(self):
        return self.bins.taxid.unique()

    def get_binids(self):
        return self.bins.binid.unique()

    def get_seqids(self):
        return self.bins.seqid.unique()

    def get_binid_length_sum(self):
        return self.bins.groupby(['binid'])['length'].sum()
        
    def get_max_bin_length(self):
        return self.get_binid_length_sum().max()

    def get_specialization_taxid(self):
        return self.bins[['specialization','taxid']].drop_duplicates().set_index('specialization')

    def merge(self, bins):
        self.bins = pd.concat([self.bins, bins.get_bins()])

    def write_acc_bin_file(self, acc_bin_file, binids: set=None):
        if binids is None:
            self.bins.to_csv(acc_bin_file, header=False, index=False, columns=['seqid','seqstart','seqend','binid'], sep='\t')
        else:
            self.bins.loc[self.bins['binid'].isin(binids)].to_csv(acc_bin_file, header=False, index=False, columns=['seqid','seqstart','seqend','binid'], sep='\t')

    def write_map_file(self, map_file, use_assembly):
        if use_assembly:
            self.bins[['specialization','binid']].drop_duplicates().to_csv(map_file,header=False, index=False, sep='\t')
        else:
            self.bins[['taxid','binid']].drop_duplicates().to_csv(map_file,header=False, index=False, sep='\t')
       
class Tax:
    def __init__(self, tax_files: list=None, ncbi_nodes: str=None, ncbi_names: str=None):
        # nodes[node] = (parent, rank, name)
        self.nodes = {}
        # default root node
        self.nodes["1"] = ("0", "no rank", "root")
        if tax_files is not None:
            self.parse(tax_files)
        elif ncbi_nodes is not None:
            self.parse_ncbi(ncbi_nodes, ncbi_names)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'Tax({})'.format(', '.join(args))

    def add_nodes(self, new_nodes, label):
        # add nodes from a pandas dataframe [group taxid]
        for node,parent in new_nodes.iterrows():
            if node not in self.nodes:
                self.nodes[node] = (parent['taxid'],label,node) # repeat node on name

    def merge(self, extra_tax): # duplicates solved by current tax
        for node in extra_tax.nodes:
            if node not in self.nodes: # if present in the current tax
                self.nodes[node] = extra_tax.nodes[node]

    def filter(self, filter_nodes):
        new_nodes = {}
        # Filter nodes for used taxids
        for node in filter_nodes:
            if node in self.nodes: # if present in the current tax
                t = node
                while t!="0": # while not at root, get lineage
                    if t in new_nodes: break # branch already added
                    new_nodes[t] = self.nodes[t] # copy node
                    t = self.nodes[t][0] # get parent
        self.nodes = new_nodes

    def parse(self, tax_files):
        for tax_file in tax_files:
            with open(tax_file , 'r') as file:
                for line in file:
                    node, parent, rank, name = line.rstrip().split("\t")
                    if node not in self.nodes:
                        self.nodes[node] = (parent,rank,name)

    def parse_ncbi(self, ncbi_nodes, ncbi_names):
        names = defaultdict(lambda: "")
        if ncbi_names is not None:
            with open(ncbi_names,'r') as file:
                for line in file:
                    node, name, _, name_class = line.split('\t|\t') # READ names -> fields (1:TAXID 2:NAME 3:UNIQUE NAME 4:NAME CLASS)
                    if name_class.replace('\t|\n','')=="scientific name":
                        names[node] = name
        with open(ncbi_nodes , 'r') as file:
            for line in file:
                node, parent, rank, _ = line.split('\t|\t', 3) # READ nodes -> fields (1:TAXID 2:PARENT_TAXID 3:RANK)
                if node not in self.nodes:
                    self.nodes[node] = (parent,rank,names[node])

    def write(self, file):
        # .tax: taxid/assembly <tab> parent taxid <tab> rank <tab> name
        tax_file = open(file,'w')
        for node,(parent, rank, name) in self.nodes.items():
            print(node, parent, rank, name, sep="\t", file=tax_file)
        tax_file.close()

    def get_rank(self, node, rank):
        t = node
        try:
            while self.nodes[t][1]!=rank and t!="1": t = self.nodes[t][0]
        except:
            return node
        return t if t!="1" else node

    def get_node_rank_fixed(self, node, fixed_ranks):
        if node!="0":
            original_rank = self.nodes[node][1]
            original_taxid = node
            while node!="0":
                if(self.nodes[node][1] in fixed_ranks):
                    #everything below species (not being assembly) is counted as species+
                    if "species+" in fixed_ranks and original_rank!="species" and original_rank!="assembly" and self.nodes[node][1]=="species":
                        return original_taxid, "species+"
                    else:
                        return node, self.nodes[node][1]
                node = self.nodes[node][0]
            return "1", "root" #no standard rank identified
        else:
            return "0", ""

if __name__ == '__main__':
    main()
