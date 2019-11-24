#!/usr/bin/env python3

import argparse, os, sys, subprocess, io, time, shlex, shutil, gzip, pickle, math, re, copy
from collections import defaultdict, OrderedDict

def main(arguments=None):

    version = '0.2.0'
    
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
    build_group_optional.add_argument('-r', '--rank',            type=str,   default='species',metavar='', help='Lowest taxonomic rank for classification [assembly,taxid,species,genus,...]. Default: species')
    build_group_optional.add_argument('-k', '--kmer-size',       type=int,   default=19,      metavar='', help='The k-mer size for the bloom filter. Default: 19')
    build_group_optional.add_argument('-n', '--hash-functions',  type=int,   default=3,       metavar='', help='The number of hash functions to use for the bloom filter. Default: 3')
    build_group_optional.add_argument('-f', '--max-fp',          type=float, default=0.05,    metavar='', help='Max. false positive rate for k-mer classification. Default: 0.05')
    build_group_optional.add_argument('-m', '--max-bloom-size',  type=int,                    metavar='', help='Approx. maximum filter size in Megabytes (MB). Will estimate best --bin-length based on --kmer-size, --hash-functions and --max-fp  [Mutually exclusive --fixed-bloom-size]')
    build_group_optional.add_argument('-l', '--bin-length',      type=int,                    metavar='', help='Maximum length (in bp) for each bin. Default: auto')
    build_group_optional.add_argument('-t', '--threads',         type=int,   default=2,       metavar='', help='Number of subprocesses/threads to use for calculations. Default: 2')
    build_group_optional.add_argument('--fixed-bloom-size',      type=int,                    metavar='', help='Fixed size for filter in Megabytes (MB), will ignore --max-fp [Mutually exclusive --max-bloom-size] ')
    build_group_optional.add_argument('--fragment-length',       type=int,   default=-1,      metavar='', help='Fragment length (in bp). Set to 0 to not fragment sequences. Default: --bin-length - --overlap-length')
    build_group_optional.add_argument('--overlap-length',        type=int,   default=300,     metavar='', help='Fragment overlap length (in bp). Should be bigger than the read length used for classification. Default: 300')
    build_group_optional.add_argument('--seq-info',              type=str, nargs="*", default=["auto"],  metavar='', help='Mode to obtain sequence information. For each sequence entry provided, ganon requires taxonomic and seq. length information. If a small number of sequences is provided (<50000) or when --rank assembly, ganon will automatically obtain data with NCBI E-utils websevices (eutils). Offline mode will download batch files from NCBI Taxonomy and look for taxonomic ids in the order provided. Options: [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot], eutils (force webservices) or auto (uses eutils or [nucl_gb nucl_wgs]). Default: auto [Mutually exclusive --seq-info-file]')
    build_group_optional.add_argument('--seq-info-file',         type=str,                               metavar='', help='Pre-generated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly id]) [Mutually exclusive --seq-info]')
    build_group_optional.add_argument('--taxdump-file',          type=str, nargs="*",                    metavar='', help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
    build_group_optional.add_argument('--input-directory',       type=str,                    metavar='', help='Directory containing input files')
    build_group_optional.add_argument('--input-extension',       type=str,                    metavar='', help='Extension of files to use with --input-directory')

    # Extra
    build_group_optional.add_argument('--verbose', default=False, action='store_true', help='Verbose mode for ganon')
    build_group_optional.add_argument('--ganon-path', type=str, default="", help=argparse.SUPPRESS)
    build_group_optional.add_argument('--taxsbp-path', type=str, default="", help=argparse.SUPPRESS)
    build_group_optional.add_argument('--n-refs', type=int, help=argparse.SUPPRESS)
    build_group_optional.add_argument('--n-batches', type=int, help=argparse.SUPPRESS)

    ####################################################################################################

    update_parser = argparse.ArgumentParser(description='Update options', add_help=False)

    # Required
    update_group_required = update_parser.add_argument_group('required arguments')
    update_group_required.add_argument('-d', '--db-prefix',         required=True,  type=str,               metavar='db_prefix',        help='Database prefix')
    update_group_required.add_argument('-i', '--input-files',       required=False, type=str, nargs="*",    metavar='',  help='Input reference sequence fasta files [.gz]')
    
    # Defaults
    update_group_optional = update_parser.add_argument_group('optional arguments')
    update_group_optional.add_argument('-o', '--output-db-prefix',                  type=str,                               metavar='', help='Alternative output database prefix. Default: overwrite current --db-prefix')
    #update_group_optional.add_argument('-c', '--update-complete',                             default=False, action='store_true', help='Update complete bins, removing sequences. Input file should be complete, not only new sequences.')
    update_group_optional.add_argument('-t', '--threads',                           type=int, default=2,                    metavar='', help='set the number of subprocesses/threads to use for calculations. Default: 2')
    update_group_optional.add_argument('--seq-info',              type=str, nargs="*", default=["auto"],  metavar='', help='Mode to obtain sequence information. For each sequence entry provided, ganon requires taxonomic and seq. length information. If a small number of sequences is provided (<50000) or when --rank assembly, ganon will automatically obtained data with NCBI E-utils websevices (eutils). Offline mode will download batch files from NCBI Taxonomy and look for taxonomic ids in the order provided. Options: [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot], eutils (force webservices) or auto (uses eutils or [nucl_gb nucl_wgs]). Default: auto [Mutually exclusive --seq-info-file]')
    update_group_optional.add_argument('--seq-info-file',         type=str,                               metavar='', help='Pre-generated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly id]) [Mutually exclusive --seq-info]')
    update_group_optional.add_argument('--taxdump-file',          type=str, nargs="*",                    metavar='', help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
    update_group_optional.add_argument('--input-directory',       type=str,                    metavar='', help='Directory containing input files')
    update_group_optional.add_argument('--input-extension',       type=str,                    metavar='', help='Extension of files to use with --input-directory')

    # Extra
    update_group_optional.add_argument('--verbose', default=False, action='store_true', help='Verbose mode for ganon')
    update_group_optional.add_argument('--ganon-path', type=str, default="", help=argparse.SUPPRESS)
    update_group_optional.add_argument('--taxsbp-path', type=str, default="", help=argparse.SUPPRESS)
    update_group_optional.add_argument('--n-refs', type=int, help=argparse.SUPPRESS)
    update_group_optional.add_argument('--n-batches', type=int, help=argparse.SUPPRESS)

    ####################################################################################################

    classify_parser = argparse.ArgumentParser(description='Classification options', add_help=False)
    
    # Required
    classify_group_required = classify_parser.add_argument_group('required arguments')
    classify_group_required.add_argument('-d', '--db-prefix',    required=True, type=str,              nargs="*", metavar='db_prefix', help='Database prefix[es]')
    classify_group_required.add_argument('-r', '--single-reads',        required=False, type=str,              nargs="*", metavar='reads.fq[.gz]', help='Multi-fastq[.gz] file[s] to classify')
    classify_group_required.add_argument('-p', '--paired-reads',        required=False, type=str,              nargs="*", metavar='reads.1.fq[.gz] reads.2.fq[.gz]', help='Multi-fastq[.gz] pairs of file[s] to classify')

    # Defaults
    classify_group_optional = classify_parser.add_argument_group('optional arguments')
    classify_group_optional.add_argument('-c', '--db-hierarchy',                type=str, default=["1"], nargs="*", metavar='int', help='Hierachy definition, one for each database input. Can also be string, but input will be always sorted (e.g. 1 1 2 3). Default: 1')
    classify_group_optional.add_argument('-m', '--min-kmers',                   type=float, default=[0.25], nargs="*", metavar='int', help='Min. percentage of k-mers matching to consider a read assigned. Can be used alternatively to --max-error for reads of variable size. Single value or one per database (e.g. 0.5 0.7 1 0.25). Default: 0.25 ')
    classify_group_optional.add_argument('-e', '--max-error',                   type=int,                nargs="*", metavar='int', help='Max. number of errors allowed. Single value or one per database (e.g. 3 3 4 0) [Mutually exclusive --min-kmers]')
    classify_group_optional.add_argument('-u', '--max-error-unique',            type=int, default=[-1],  nargs="*", metavar='int', help='Max. number of errors allowed for unique assignments after filtering. Matches below this error rate will not be discarded, but assigned to parent taxonomic level. Single value or one per hierachy (e.g. 0 1 2). -1 to disable. Default: -1')    
    classify_group_optional.add_argument('-f', '--offset',                      type=int, default=2,              metavar='', help='Number of k-mers to skip during clasification. Can speed up analysis but may reduce recall. (e.g. 1 = all k-mers, 3 = every 3rd k-mer). Function must be enabled on compilation time with -DGANON_OFFSET=ON. Default: 2')
    classify_group_optional.add_argument('-o', '--output-file-prefix',          type=str, default="",             metavar='', help='Output file name prefix: .out for complete results / .lca for LCA results / .rep for report. Empty to print to STDOUT (only with lca). Default: ""')
    classify_group_optional.add_argument('-n', '--output-unclassified-file',    type=str, default="",             metavar='', help='Output file for unclassified reads headers. Empty to not output. Default: ""')
    classify_group_optional.add_argument('-s', '--split-output-file-hierarchy', default=False, action='store_true',               help='Split output in multiple files by hierarchy. Appends "_hierachy" to the --output-file definiton.')
    classify_group_optional.add_argument('--skip-lca',                    default=False, action='store_true',               help='Skip LCA step and output multiple matches. --max-error-unique will not be applied')
    classify_group_optional.add_argument('--skip-reports',                default=False, action='store_true',               help='Skip reports')
    classify_group_optional.add_argument('-k', '--ranks',                       type=str, default=[],   nargs="*",                help='Ranks for the final report. "all" for all indentified ranks. empty for default ranks: superkingdom phylum class order family genus species species+ assembly')
    classify_group_optional.add_argument('-t', '--threads',                     type=int, default=3,              metavar='', help='Number of subprocesses/threads. Default: 3)')
    # Extra
    classify_group_optional.add_argument('--verbose',                           default=False, action='store_true',  help='Output in verbose mode for ganon-classify')
    classify_group_optional.add_argument('--ganon-path', type=str, default="", help=argparse.SUPPRESS)
    classify_group_optional.add_argument('--n-reads', type=int, help=argparse.SUPPRESS)
    classify_group_optional.add_argument('--n-batches', type=int, help=argparse.SUPPRESS)

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

    args = parser.parse_args()

    if len(sys.argv[1:])==0: # Print help calling script without parameters
        parser.print_help() 
        return 0
    
    tx_total = time.time()

    args.ganon_path = args.ganon_path + "/" if args.ganon_path else ""

    if args.which=='build' or args.which=='update':
        # if path is given, look for binaries only there
        ganon_build_paths = [args.ganon_path, args.ganon_path+"build/"] if args.ganon_path else [None, "build/"]
        for p in ganon_build_paths:
            ganon_build_exec = shutil.which("ganon-build", path=p)
            if ganon_build_exec is not None: break
        if ganon_build_exec is None:
            print_log("ganon-build binary was not found. Please inform a specific path with --ganon-path\n")
            return 1

        ganon_get_len_taxid_paths = [args.ganon_path, args.ganon_path+"scripts/", args.ganon_path+"../scripts/"] if args.ganon_path else [None, "scripts/"]
        for p in ganon_get_len_taxid_paths:
            ganon_get_len_taxid_exec = shutil.which("ganon-get-len-taxid.sh", path=p)
            if ganon_get_len_taxid_exec is not None: break
        if ganon_get_len_taxid_exec is None:
            print_log("ganon-get-len-taxid.sh script was not found. Please inform a specific path with --ganon-path\n")
            return 1

        taxsbp_paths = [args.taxsbp_path + "/"] if args.taxsbp_path else [None, "taxsbp/"]
        for p in taxsbp_paths:
            taxsbp_exec = shutil.which("taxsbp", path=p)
            if taxsbp_exec is not None: break
            taxsbp_exec = shutil.which("taxsbp.py", path=p)
            if taxsbp_exec is not None: break
        if taxsbp_exec is None:
            print_log("TaxSBP script (taxsbp or taxsbp.py) were not found. Please inform the path of the scripts with --taxsbp-path\n")
            return 1

        if args.taxdump_file and ((len(args.taxdump_file)==1 and not args.taxdump_file[0].endswith(".tar.gz")) or len(args.taxdump_file)>3):
            print_log("Please provide --taxdump-file taxdump.tar.gz or --taxdump-file nodes.dmp names.dmp [merged.dmp] or leave it empty for automatic download \n")
            return 1

        if not args.input_files and not args.input_directory:
            print_log("Please provide files with --input-files and/or --input-directory with --input-extension \n")
            return 1
        elif args.input_directory and not args.input_extension:
            print_log("Please provide the --input-extension when using --input-directory \n")
            return 1
        elif args.input_directory and not args.input_files:
            args.input_files = [] # initializate for adding files later

        db_prefix = args.db_prefix
        output_folder = os.path.abspath(os.path.dirname(db_prefix)) + "/"
        tmp_output_folder = db_prefix + "_tmp/"
        db_prefix_ibf = db_prefix + ".ibf"
        db_prefix_map = db_prefix + ".map"
        db_prefix_tax = db_prefix + ".tax"
        db_prefix_gnn = db_prefix + ".gnn"

    if args.which=='build': #If set (!=0), should be smaller than fragment
        if args.fragment_length>0 and args.overlap_length > args.fragment_length:
            print_log("--overlap-length cannot be bigger than --fragment-length\n")
            return 1

        if args.fixed_bloom_size and not args.bin_length:
            print_log("please set the --bin-length to use --fixed-bloom-size\n")
            return 1

    if args.which=='classify':

        ganon_classify_paths = [args.ganon_path, args.ganon_path+"build/"] if args.ganon_path else [None, "build/"]
        for p in ganon_classify_paths:
            ganon_classify_exec = shutil.which("ganon-classify", path=p)
            if ganon_classify_exec is not None: break
        if ganon_classify_exec is None:
            print_log("ganon-classify binary was not found. Please inform a specific path with --ganon-path\n")
            return 1
        if not args.single_reads and not args.paired_reads:
            print_log("Please provide reads with the paramenter -r/--single-reads and/or -p/--paired-reads\n")
            return 1

        if args.skip_lca and not args.output_file_prefix:
            print_log("--output-file-prefix is mandatory without LCA\n")
            return 1

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    if args.which=='build':     
        use_assembly=True if args.rank=="assembly" else False
        taxsbp_input_file, ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = prepare_files(args, tmp_output_folder, use_assembly, ganon_get_len_taxid_exec)

        # Set bin length
        if args.bin_length: # user defined
            bin_length = args.bin_length
        else:
            tx = time.time()
            print_log("Estimating best bin length... ")
            bin_length = estimate_bin_len(args, taxsbp_input_file, ncbi_nodes_file, use_assembly)
            if bin_length==0: return 1
            print_log(str(bin_length) + "bp. ")
            print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

        # Set fragment length
        if args.fragment_length==-1: # if ==-1 set default
            fragment_length = bin_length - args.overlap_length
        elif args.fragment_length==0: # if ==0 deactivate
            fragment_length = 0
        else: # user input
            fragment_length = args.fragment_length - args.overlap_length 

        tx = time.time()
        print_log("Running taxonomic clustering (TaxSBP)... ")
        run_taxsbp_cmd = " ".join([taxsbp_exec,
                                   "-f " + taxsbp_input_file,
                                   "-n " + ncbi_nodes_file,
                                   "-m " + ncbi_merged_file if ncbi_merged_file else "",
                                   "-l " + str(bin_length),
                                   "-r " + ("assembly" if use_assembly else args.rank),
                                   "-z " + "assembly" if use_assembly else "",
                                   "-a " + str(fragment_length) if fragment_length else "",
                                   "-o " + str(args.overlap_length) if fragment_length else ""])
        stdout, stderr, errcode = run(run_taxsbp_cmd, print_stderr=True)
        
        acc_bin_file = tmp_output_folder + "acc_bin.txt"
        bins, actual_number_of_bins, unique_taxids, max_length_bin, group_taxid, _ = taxsbp_output_files(stdout, acc_bin_file, db_prefix_map, use_assembly, fragment_length)
        print_log(str(actual_number_of_bins) + " bins created. ")
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

        # aproximate number of unique k-mers by just considering that they are all unique
        max_kmer_count = max_length_bin-args.kmer_size+1
        print_log("Approximate (upper-bound) # unique k-mers: " + str(max_kmer_count) + "\n")
        # get optimal number of bins to get correct fp rate based on IBF implementation (always next multiple of 64)
        optimal_number_of_bins = optimal_bins(actual_number_of_bins)

        # define bloom filter size based on given false positive
        MBinBits = 8388608
        if not args.fixed_bloom_size:
            bin_size_bits = math.ceil(-(1/((1-args.max_fp**(1/float(args.hash_functions)))**(1/float(args.hash_functions*max_kmer_count))-1)))   
            print_log("Bloom filter calculated size with fp<=" + str(args.max_fp) + ": " + str("{0:.4f}".format((bin_size_bits*optimal_number_of_bins)/MBinBits)) + "MB (" + str(bin_size_bits) + " bits/bin * " + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])\n")
        else:
            bin_size_bits = math.ceil((args.fixed_bloom_size * MBinBits)/optimal_number_of_bins);
            estimated_max_fp = (1-((1-(1/float(bin_size_bits)))**(args.hash_functions*max_kmer_count)))**args.hash_functions
            print_log("Bloom filter calculated max. fp with size=" + str(args.fixed_bloom_size) + "MB: " + str("{0:.4f}".format(estimated_max_fp) + " ("  + str(optimal_number_of_bins) + " optimal bins [" + str(actual_number_of_bins) + " real bins])\n"))

        tx = time.time()
        print_log("Building index (ganon-build)... \n")
        run_ganon_build_cmd = " ".join([ganon_build_exec,
                                        "--seqid-bin-file " + acc_bin_file,
                                        "--filter-size-bits " + str(bin_size_bits*optimal_number_of_bins) if args.max_fp else "--filter-size " + str(args.fixed_bloom_size),
                                        "--kmer-size " + str(args.kmer_size),
                                        "--hash-functions " + str(args.hash_functions),
                                        "--threads " + str(args.threads),
                                        "--output-filter-file " + db_prefix_ibf,
                                        "--verbose" if args.verbose else "",
                                        "--n-refs " + str(args.n_refs) if args.n_refs is not None else "",
                                        "--n-batches " + str(args.n_batches) if args.n_batches is not None else "",
                                        "--reference-files " + ",".join([file for file in args.input_files]) if args.input_files else ""
                                        "--directory-reference-files " + args.input_directory if args.input_directory else ""
                                        "--extension " + args.input_extension if args.input_extension else ""])
        stdout, stderr, errcode = run(run_ganon_build_cmd, print_stderr=True)
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

        tx = time.time()
        print_log("Building database files... ")

        # Build .tax file
        build_tax(db_prefix_tax, unique_taxids, group_taxid, ncbi_nodes_file, ncbi_names_file, use_assembly)

        gnn = Gnn(kmer_size=args.kmer_size, 
                hash_functions=args.hash_functions, 
                number_of_bins=actual_number_of_bins, 
                rank=args.rank,
                bin_length=bin_length,
                fragment_length=fragment_length,
                overlap_length=args.overlap_length,
                bins=bins)
        gnn.write(db_prefix_gnn)

        # Delete temp files
        shutil.rmtree(tmp_output_folder)

        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
        
 
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    elif args.which=='update':  
        tx = time.time()

        tmp_db_prefix = tmp_output_folder + "tmp"
        tmp_db_prefix_ibf = tmp_db_prefix + ".ibf"
        tmp_db_prefix_gnn = tmp_db_prefix + ".gnn"
        tmp_db_prefix_map = tmp_db_prefix + ".map"
        tmp_db_prefix_tax = tmp_db_prefix + ".tax"

        gnn = Gnn(file=db_prefix_gnn)

        use_assembly=True if gnn.rank=="assembly" else False

        taxsbp_input_file, ncbi_nodes_file, ncbi_merged_file, ncbi_names_file = prepare_files(args, tmp_output_folder, use_assembly, ganon_get_len_taxid_exec)
        
        bins_file = tmp_output_folder + "ganon.bins"
        gnn.write_bins(bins_file)

        print_log("Running taxonomic clustering (TaxSBP)... \n")
        run_taxsbp_cmd = " ".join([taxsbp_exec,
                                   "-u " + bins_file,
                                   "-f " + taxsbp_input_file,
                                   "-n " + ncbi_nodes_file,
                                   "-m " + ncbi_merged_file if ncbi_merged_file else "",
                                   "-l " + str(gnn.bin_length),
                                   "-r " + ("assembly" if use_assembly else gnn.rank),
                                   "-z " + "assembly" if use_assembly else "",
                                   "-a " + str(gnn.fragment_length) if gnn.fragment_length else "",
                                   "-o " + str(gnn.overlap_length) if gnn.fragment_length else ""])
        stdout, stderr, errcode = run(run_taxsbp_cmd, print_stderr=True)

        acc_bin_file = tmp_output_folder + "acc_bin.txt"
        bins, number_of_updated_bins, unique_taxids, max_length_bin, group_taxid, last_bin = taxsbp_output_files(stdout, acc_bin_file, tmp_db_prefix_map, use_assembly, gnn.fragment_length)
        number_of_new_bins = int(last_bin) + 1 - gnn.number_of_bins
        print_log(str(number_of_updated_bins) + " bins updated, " + str(number_of_new_bins) + " new. ")
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

        tx = time.time()
        print_log("Updating index (ganon-build)... \n")
        run_ganon_build_cmd = " ".join([ganon_build_exec,
                                        "--update-filter-file " + db_prefix_ibf,
                                        "--seqid-bin-file " + acc_bin_file,
                                        "--output-filter-file " + tmp_db_prefix_ibf,
                                        "--threads " + str(args.threads),                                
                                        "--verbose" if args.verbose else "",
                                        "--n-refs " + str(args.n_refs) if args.n_refs is not None else "",
                                        "--n-batches " + str(args.n_batches) if args.n_batches is not None else "",
                                        "--reference-files " + ",".join([file for file in args.input_files]) if args.input_files else ""
                                        "--directory-reference-files " + args.input_directory if args.input_directory else ""
                                        "--extension " + args.input_extension if args.input_extension else ""])

        stdout, stderr, errcode = run(run_ganon_build_cmd, print_stderr=True)
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

        tx = time.time()
        print_log("Updating database files ... ")

        # Build .tax file
        build_tax(tmp_db_prefix_tax, unique_taxids, group_taxid, ncbi_nodes_file, ncbi_names_file, use_assembly)

        # move IBF
        shutil.move(tmp_db_prefix_ibf, args.output_db_prefix + ".ibf" if args.output_db_prefix else db_prefix_ibf)
        
        # write GNN in a different location
        gnn.number_of_bins+=number_of_new_bins # add new bins count
        gnn.bins.extend(bins) # save new bins from taxsbp
        gnn.write(args.output_db_prefix + ".gnn" if args.output_db_prefix else db_prefix_gnn)

        # update and write TAX
        update_tax(db_prefix_tax, tmp_db_prefix_tax, args.output_db_prefix + ".tax" if args.output_db_prefix else db_prefix_tax)
        
        # update and write MAP
        update_map(db_prefix_map, tmp_db_prefix_map, args.output_db_prefix + ".map" if args.output_db_prefix else db_prefix_map)

        # Delete temp files
        shutil.rmtree(tmp_output_folder)

        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

    elif args.which=='classify':

        if not args.output_file_prefix:
            ganon_classify_output_file = "ganon-classify-tmp.out" 
            args.skip_reports=True #output only stdout, no reports
        else:
            ganon_classify_output_file = args.output_file_prefix+".out"
        
        output_hierarchy = {}
        # if there's LCA and output files generate hierarchy
        if not args.skip_lca:
            if not args.output_file_prefix:
                output_hierarchy[None] = {'db_prefixes': args.db_prefix, 'out_file': ganon_classify_output_file, 'lca_file': "", 'rep_file': ""}
            else:
                # build hierarchy structure for output files
                if len(args.db_hierarchy) > 1 and args.split_output_file_hierarchy:
                    for dbid,dbp in enumerate(args.db_prefix):
                        hierarchy_name = args.db_hierarchy[dbid]
                        h_prefix = "_"+hierarchy_name
                        if hierarchy_name not in output_hierarchy: 
                            output_hierarchy[hierarchy_name] = {'db_prefixes': [], 'out_file': args.output_file_prefix+".out"+h_prefix, 'lca_file': args.output_file_prefix+".lca"+h_prefix, 'rep_file': args.output_file_prefix+".rep"+h_prefix}
                            if os.path.exists(output_hierarchy[hierarchy_name]['out_file']): os.remove(output_hierarchy[hierarchy_name]['out_file'])
                            if os.path.exists(output_hierarchy[hierarchy_name]['lca_file']): os.remove(output_hierarchy[hierarchy_name]['lca_file'])
                            if os.path.exists(output_hierarchy[hierarchy_name]['rep_file']): os.remove(output_hierarchy[hierarchy_name]['rep_file'])
                        output_hierarchy[hierarchy_name]['db_prefixes'].append(dbp)
                else: # no split or no hierarchy, output together
                    output_hierarchy[None] = {'db_prefixes': args.db_prefix, 'out_file': args.output_file_prefix+".out", 'lca_file': args.output_file_prefix+".lca", 'rep_file': args.output_file_prefix+".rep"}
                
                # sort it to output in the same order as ganon-classify
                output_hierarchy = OrderedDict(sorted(output_hierarchy.items(), key=lambda t: t[0]))
        
        tx = time.time()
        print_log("Classifying reads (ganon-classify)... \n")
        run_ganon_classify = " ".join([ganon_classify_exec,
                                       "--bloom-filter " + ",".join([db_prefix+".filter" for db_prefix in args.db_prefix]),
                                       "--group-bin " + ",".join([db_prefix+".map" for db_prefix in args.db_prefix]),
                                       "--filter-hierarchy " + ",".join([str(h) for h in args.db_hierarchy]),
                                       "--max-error " + ",".join([str(me) for me in args.max_error]) if not args.min_kmers else "--min-kmers " + ",".join([str(mk) for mk in args.min_kmers]),
                                       "--max-error-unique " + ",".join([str(meu) for meu in args.max_error_unique]),
                                       "--threads " + str(args.threads),
                                       "--offset " + str(args.offset),
                                       "--output-file " + ganon_classify_output_file,
                                       "--split-output-file-hierarchy " if args.split_output_file_hierarchy else "",
                                       "--output-unclassified-file " + args.output_unclassified_file if args.output_unclassified_file else "",
                                       "--verbose" if args.verbose else "",
                                       "--n-reads " + str(args.n_reads) if args.n_reads is not None else "",
                                       "--n-batches " + str(args.n_batches) if args.n_batches is not None else "",
                                       "--single-reads " +  ",".join(args.single_reads) if args.single_reads else "",
                                       "--paired-reads " +  ",".join(args.paired_reads) if args.paired_reads else ""])
        process_ganon_classify = run(run_ganon_classify, blocking=False)
        
        if not args.skip_lca:
            try:
                from pylca.pylca import LCA
            except ImportError:
                print_log("LCA module not found (pylca)")
                args.skip_lca=True

        if not args.skip_lca:
            reports = {}
            merged_filtered_nodes = {}

            for hierarchy_name, hierarchy in output_hierarchy.items():

                # Merge group-taxid information and nodes from multiple databases
                merged_filtered_nodes[hierarchy_name] = {}
                for db_prefix in hierarchy["db_prefixes"]:
                    _, _, _, rank, _, _, _, filtered_nodes = parse_db_prefix_gnn(db_prefix+".gnn")
                    use_assembly=True if rank=="assembly" else False # if one of them uses assembly should be True
                    merged_filtered_nodes[hierarchy_name].update(filtered_nodes)
                
                # pre build LCA with used nodes
                L = LCA({tx:parent for tx,(parent,_,_) in merged_filtered_nodes[hierarchy_name].items()})

                # wait ganon-classify to create the file
                while not os.path.isfile(hierarchy["out_file"]): time.sleep(1)

                # redirect output for lca
                if args.output_file_prefix: sys.stdout = open(hierarchy["lca_file"],'w')

                reports[hierarchy_name] = generate_lca_rep(args, hierarchy["out_file"], merged_filtered_nodes[hierarchy_name], L, use_assembly)

                # Rm file if tmp
                if not args.output_file_prefix: 
                    os.remove(ganon_classify_output_file) 
                else: # close open file if not
                    sys.stdout.close()
                    sys.stdout = sys.__stdout__ #return to stdout
                
        stdout, stderr = process_ganon_classify.communicate()
        print_log(stderr)
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

        if not args.skip_lca and not args.skip_reports:
            tx = time.time()
            print_log("Generating reports... ")

            # get reads classified and unclassified from strerr
            re_out = re.search(r"\d+\ssequences classified", stderr)
            seq_cla = int(re_out.group().split(" ")[0]) if re_out is not None else 0
            re_out = re.search(r"\d+\ssequences unclassified", stderr)
            seq_unc = int(re_out.group().split(" ")[0]) if re_out is not None else 0
            total_reads = seq_cla+seq_unc

            for hierarchy_name, hierarchy in output_hierarchy.items():
                
                with open(hierarchy["rep_file"], 'w') as rfile:
                    rfile.write("unclassified" +"\t"+ str(seq_unc) +"\t"+ str("%.5f" % ((seq_unc/total_reads)*100)) +"\t"+ "0" +"\t"+ "0" +"\t"+ "0" +"\t"+ "-" +"\t"+ "-" + "\n")
                    for assignment in sorted(reports[hierarchy_name], key=lambda k: reports[hierarchy_name][k]['count'], reverse=True):
                        rfile.write(assignment +"\t"+ str(reports[hierarchy_name][assignment]['count']) +"\t"+ str("%.5f" % ((reports[hierarchy_name][assignment]['count']/total_reads)*100)) +"\t"+ str(reports[hierarchy_name][assignment]['assignments']) +"\t"+ str(reports[hierarchy_name][assignment]['unique']) +"\t"+ merged_filtered_nodes[hierarchy_name][assignment][2] +"\t"+ merged_filtered_nodes[hierarchy_name][assignment][1] + "\n")

            final_report_file = args.output_file_prefix+".tre"
            print_final_report(reports, merged_filtered_nodes, seq_unc, total_reads, final_report_file, args.ranks)
            print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
        
    print_log("Total elapsed time: " + str("%.2f" % (time.time() - tx_total)) + " seconds.\n")


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def generate_lca_rep(args, out_file, merged_filtered_nodes_hierarchy, L, use_assembly):
    rep = defaultdict(lambda: {'count': 0, 'assignments': 0, 'unique': 0})
    with open(out_file, "r") as file:
        # read first line
        old_readid, cl, kc, done = get_output_line(file)
        if not done:
            assignments = set([cl])
            max_kmer_count = int(kc)

        while not done:
            readid, cl, kc, done = get_output_line(file)
            # next read matches, print old
            # if done=True, readid=="" to print last entries
            if readid != old_readid: 
                len_assignments = len(assignments)
                taxid_lca = get_lca_read(copy.copy(assignments), max_kmer_count, merged_filtered_nodes_hierarchy, L, use_assembly)
                print(old_readid, taxid_lca, max_kmer_count, sep="\t")
                if not args.skip_reports:
                    rep[taxid_lca]['count'] += 1
                    rep[taxid_lca]['assignments'] += len_assignments
                    if len_assignments==1 and max_kmer_count>0: rep[taxid_lca]['unique'] += 1 # only count as unique if was not negative (meaning it didn't pass unique error filter)
                        
                assignments.clear()
                max_kmer_count=0
        
            if not done: #if not last line
                assignments.add(cl)
                kc = int(kc)
                # account for unique filtering with abs but keep it negative to retain information on output
                max_kmer_count = kc if abs(kc) > abs(max_kmer_count) else max_kmer_count 
                old_readid = readid 
    return rep


def get_lca_read(assignments, max_kmer_count, merged_filtered_nodes_hierarchy, L, use_assembly):
    lca_taxid=""
    if len(assignments)==1: # unique match or same taxid (but one or more assemblies)  
        if max_kmer_count<0: # unique matches in ganon-classify, get leaf taxid (assembly) or parent taxid
            lca_taxid = merged_filtered_nodes_hierarchy[assignments.pop()][0]
        else:
            lca_taxid = assignments.pop()
    else:
        if use_assembly: # get taxids from assembly (reduce number of entries, could be done on the LCA, but is done here for speed-up)
            assignments = set(map(lambda x: merged_filtered_nodes_hierarchy[x][0], assignments))  # Recover taxids from assignments (assembly)
            if len(assignments)==1: # If all assignments are on the same taxid, no need to lca and return it
                lca_taxid = assignments.pop()
    
    if lca_taxid:
        return lca_taxid
    else:
        # Instead of return the LCA of all pairs L(*assignments), do it in pairs and saving in a set (no redudancy)
        # reducing the number of calculations by 2 in the best case (around 25% less in some tests)
        lca_taxid = set()
        while(True):
            if len(assignments)<=1: 
                if len(assignments)==1: lca_taxid.add(assignments.pop()) #odd number of assignments
                if len(lca_taxid)==1: 
                    break
                else:
                    assignments = lca_taxid
                    lca_taxid = set()
            lca_taxid.add(L(assignments.pop(),assignments.pop()))
        return lca_taxid.pop()

def prepare_files(args, tmp_output_folder, use_assembly, ganon_get_len_taxid_exec):
    # Create temporary working directory
    if os.path.exists(tmp_output_folder): shutil.rmtree(tmp_output_folder) # delete if already exists
    os.makedirs(tmp_output_folder)
    
    # get files from
    if args.input_directory and args.input_extension:
        for file in os.listdir(args.input_directory):
            if file.endswith(args.input_extension):
                args.input_files.append(os.path.join(args.input_directory, file))

    # Prepare TaxSBP input file
    if args.seq_info_file: # file already provided 
        taxsbp_input_file = args.seq_info_file
    else: # retrieve info
        taxsbp_input_file = retrieve_ncbi(tmp_output_folder, args.input_files, args.threads, ganon_get_len_taxid_exec, args.seq_info, use_assembly)

    if not args.taxdump_file:
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(get_taxdump(tmp_output_folder), tmp_output_folder)
    elif args.taxdump_file[0].endswith(".tar.gz"):
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(args.taxdump_file[0], tmp_output_folder)
    else:
        ncbi_nodes_file = args.taxdump_file[0]
        ncbi_names_file = args.taxdump_file[1]
        ncbi_merged_file =  args.taxdump_file[2] if len(args.taxdump_file)==3 else ""

    return taxsbp_input_file, ncbi_nodes_file, ncbi_merged_file, ncbi_names_file

def run(cmd, output_file=None, print_stderr=False, shell=False, blocking=True):
    try:
        errcode=0
        process = subprocess.Popen(shlex.split(cmd) if not shell else cmd, shell=shell, universal_newlines=True, stdout=subprocess.PIPE if output_file is None else open(output_file, 'w'), stderr=subprocess.PIPE)   
        if blocking:
            stdout, stderr = process.communicate() # wait for the process to terminate
            errcode = process.returncode
            if errcode!=0: raise Exception()
            if print_stderr: print_log(stderr)
            return stdout, stderr, errcode
        else:
            return process
    #except OSError as e: # The most common exception raised is OSError. This occurs, for example, when trying to execute a non-existent file. Applications should prepare for OSError exceptions.
    #except ValueError as e: #A ValueError will be raised if Popen is called with invalid arguments.
    except Exception as e:
        print_log('The following command failed to execute:\n'+cmd)
        print_log(str(e)+"\n")
        print_log("Errorcode: "+str(errcode)+"\n")
        print_log("Error: "+stderr+"\n")
        raise

def taxsbp_output_files(stdout, acc_bin_file, db_prefix_map, use_assembly, fragment_length):
    bins = []
    bin_group = dict() # for unique map entries
    unique_taxids = set() # for nodes
    bin_length = defaultdict(int) # calculate max bin length
    group_taxid = {} # set group taxid connection for nodes in case of assembly

    acc_bin = open(acc_bin_file,'w') # input for ganon-build for build
    for line in stdout.split("\n"):
        if line:
            #acc, length, taxid, [group/rank taxid,] binno
            if use_assembly:
                acc, length, taxid, group, binno = line.split("\t")
                group_taxid[group] = taxid
            else:
                acc, length, taxid, binno = line.split("\t")
                group = taxid # taxid is the classification group

            # if sequences are fragmentes
            if fragment_length: 
                acc, coord = acc.split("/")
                frag_start, frag_end = coord.split(":")
            else:
                frag_start = 1
                frag_end = length

            bin_length[int(binno)]+=int(length)
            bin_group[binno] = group # save unique mapping
            unique_taxids.add(taxid)
            bins.append(line)
            print(acc, frag_start, frag_end, binno, sep="\t", file=acc_bin)
    acc_bin.close()

    # write .map
    with open(db_prefix_map, 'w') as file:
        for binno, group in bin_group.items(): 
            file.write(group + "\t" + binno + "\n")

    return bins, len(bin_length), unique_taxids, max(bin_length.values()), group_taxid, max(bin_length.keys())

def get_taxdump(tmp_output_folder):
    tx = time.time()
    print_log("Downloading taxdump... ")
    taxdump_file = tmp_output_folder+'taxdump.tar.gz'
    run_wget_taxdump_cmd = 'wget -qO {0} "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"'.format(taxdump_file)
    stdout, stderr, errcode = run(run_wget_taxdump_cmd, print_stderr=True)
    print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
    return taxdump_file

def unpack_taxdump(taxdump_file, tmp_output_folder):
    unpack_taxdump_cmd = 'tar xf {0} -C "{1}" nodes.dmp merged.dmp names.dmp'.format(taxdump_file, tmp_output_folder)
    stdout, stderr, errcode = run(unpack_taxdump_cmd, print_stderr=True)
    return tmp_output_folder+'nodes.dmp', tmp_output_folder+'names.dmp', tmp_output_folder+'merged.dmp'

def get_accession2taxid(acc2txid, tmp_output_folder):
    tx = time.time()
    acc2txid_file = acc2txid + ".accession2taxid.gz"
    print_log("Downloading " + acc2txid_file + "... ")
    acc2txid_file = tmp_output_folder + acc2txid_file
    run_wget_acc2txid_file_cmd = 'wget -qO {0} "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/{1}.accession2taxid.gz"'.format(acc2txid_file, acc2txid)
    stdout, stderr, errcode = run(run_wget_acc2txid_file_cmd, print_stderr=True)
    print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
    return acc2txid_file

def retrieve_ncbi(tmp_output_folder, files, threads, ganon_get_len_taxid_exec, seq_info, use_assembly):

    taxsbp_input_file = tmp_output_folder + 'acc_len_taxid.txt'

    # force eutils for assembly mode
    if use_assembly: seq_info = ["eutils"]

    # default accession2taxid files
    default_acc2txid = ["nucl_gb", "nucl_wgs"]

    if seq_info[0] in ["auto","eutils"]: # if using offline mode with acc2txid, get accession when extracting lenghts
        tx = time.time()
        accessions_file = tmp_output_folder + 'accessions.txt'
        print_log("Extracting accessions... ")
        # cat | zcat | gawk -> compability with osx
        run_get_header = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{print substr($1,2)}}'".format(" ".join(files), "| zcat" if files[0].endswith(".gz") else "")
        stdout, stderr, errcode = run(run_get_header, print_stderr=False, shell=True)
        count = stdout.count('\n')
        print_log(str(count) + " accessions retrieved. ")
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
        if seq_info[0]=="auto" and count>50000: 
            seq_info = default_acc2txid
        else:
            seq_info = ["eutils"]
            with open(accessions_file, "w") as af: # write accessions in file
                af.write(stdout)
        del stdout
    
    if seq_info[0]=="eutils":
        tx = time.time()
        print_log("Retrieving sequence lengths and taxid from NCBI E-utils... ")
        run_get_len_taxid_cmd = '{0} -i {1} {2}'.format(
                                ganon_get_len_taxid_exec,
                                accessions_file,
                                "-a" if use_assembly else "")
        stdout, stderr, errcode = run(run_get_len_taxid_cmd, output_file=taxsbp_input_file, print_stderr=True)
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
    
    else:
        import pandas as pd

        acc2txid_options = ["nucl_gb","nucl_wgs","nucl_est","nucl_gss","pdb","prot","dead_nucl","dead_wgs","dead_prot"]

        tx = time.time()
        print_log("Extracting sequence lengths... ")
        accessions_lengths = {}
        acc_count = 0
        for file in files:
            # cat | zcat | gawk -> compability with osx
            run_get_length = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{if (seqlen){{print seqlen}}; printf substr($1,2)\"\\t\";seqlen=0;next;}} {{seqlen+=length($0)}}END{{print seqlen}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr, errcode = run(run_get_length, print_stderr=False, shell=True)
            for line in stdout.rstrip().split("\n"):
                a, l = line.split("\t")
                if l != "0": accessions_lengths[a] = l
                acc_count+=1
        print_log(str(len(accessions_lengths)) + " entries retrieved. ")
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
        del stdout

        if len(accessions_lengths) < acc_count:
            print_log("Could not retrieve lenght for " + str(acc_count - len(accessions_lengths)) + " accessions\n")

        acc_count = len(accessions_lengths) # new count only with valid ones

        accessions_taxids = pd.DataFrame(columns=['acc','taxid'])
        for acc2txid in seq_info:
            if acc2txid not in acc2txid_options:
                print_log(acc2txid +  " is not a valid option \n")
                break
            acc2txid_file = get_accession2taxid(acc2txid, tmp_output_folder)
            tx = time.time()
            print_log("Extracting taxids from " + acc2txid_file.split("/")[-1] + "... ")
            tmp_accessions_taxids = pd.read_csv(acc2txid_file, sep='\t', header=None, skiprows=1, usecols=[1,2], names=['acc','taxid'], converters={'acc':lambda x: x if x in accessions_lengths.keys() else ""}, dtype={'taxid': str})
            tmp_accessions_taxids = tmp_accessions_taxids[tmp_accessions_taxids['acc']!=""] #keep only accessions used
            tmp_accessions_taxids = tmp_accessions_taxids[tmp_accessions_taxids['taxid']!="0"] # filter out taxid==0
            print_log(str(tmp_accessions_taxids.shape[0]) + " entries found. ")
            accessions_taxids = pd.concat([accessions_taxids,tmp_accessions_taxids])
            print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")
            if accessions_taxids.shape[0] == acc_count: #if already found all taxid for hte accessions (no need to parse all files)
                break
        
        del tmp_accessions_taxids

        if accessions_taxids.shape[0] < acc_count:
            print_log("Could not retrieve taxid for " + str(acc_count - accessions_taxids.shape[0]) + " accessions\n")

        # write file
        tx = time.time()
        print_log("Writing acc_len_taxid file... ")
        taxsbp_input = open(taxsbp_input_file,"w")
        for index, row in accessions_taxids.iterrows():
            taxsbp_input.write(row["acc"] + "\t" + accessions_lengths[row["acc"]] + "\t" + row["taxid"] + "\n")
        taxsbp_input.close()
        print_log("Done. Elapsed time: " + str("%.2f" % (time.time() - tx)) + " seconds.\n")

    return taxsbp_input_file

def update_map(old_map, new_map, output_map):
    bin_group = dict() # for unique map entries
    with open(old_map,'r') as file:
        for line in file:
            group, binno = line.rstrip().split('\t')
            bin_group[binno] = group # save unique mapping
    with open(new_map,'r') as file:
        for line in file:
            group, binno = line.rstrip().split('\t')
            bin_group[binno] = group# save unique mapping
    # write final .map
    with open(output_map, 'w') as file:
        for binno, group in bin_group.items(): 
            file.write(group + "\t" + binno + "\n")

def update_tax(old_tax, new_tax, output_tax):
    printed_nodes = set()
    tax = open(output_tax,'w')

    # copy new tax to the final output, saving nodes
    with open(new_tax,'r') as file:
        for line in file:
            taxid, _ = line.split('\t', 1)
            printed_nodes.add(taxid)
            tax.write(line)

    # merge old taxa nodes if they were not yet present - duplicates solved by new nodes
    with open(old_tax,'r') as file:
        for line in file:
            taxid, _ = line.split('\t', 1)
            if taxid not in printed_nodes:
                tax.write(line)

    tax.close()

def build_tax(db_prefix_tax, unique_taxids, group_taxid, ncbi_nodes_file, ncbi_names_file, use_assembly):

    # parse nodes.dmp and names.dmp
    ncbi_nodes, ncbi_ranks = read_nodes(ncbi_nodes_file)
    ncbi_names = read_names(ncbi_names_file)

    # .tax: taxid/assembly <tab> parent taxid <tab> rank <tab> name
    tax = open(db_prefix_tax,'w')

    # keep track of printed nodes to not repeat lineage
    printed_nodes = set()

    # if using assembly, add group to nodes - no name, repeat group
    if use_assembly:
        for group,taxid in group_taxid.items():
            print(group, taxid, "assembly", group, sep="\t", file=tax)
            printed_nodes.add(group)

    # Filter nodes for used taxids
    for leaf_taxid in unique_taxids:
        t = leaf_taxid
        while t!="0":
            if t in printed_nodes: break # branch already in the dict
            printed_nodes.add(t)
            print(t, ncbi_nodes[t], ncbi_ranks[t], ncbi_names[t], sep="\t", file=tax)
            t = ncbi_nodes[t]

    tax.close()

def read_nodes(nodes_file):
    # READ nodes -> fields (1:TAXID 2:PARENT_TAXID 3:RANK)
    nodes = {}
    ranks = {}
    with open(nodes_file,'r') as fnodes:
        for line in fnodes:
            taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
            nodes[taxid] = parent_taxid
            ranks[taxid] = rank
    nodes["1"] = "0"
    return nodes, ranks

def read_names(names_file):
    # READ names -> fields (1:TAXID 2:NAME 3:UNIQUE NAME 4:NAME CLASS)
    names = {}
    with open(names_file,'r') as fnames:
        for line in fnames:
            taxid, name, _, name_class = line.split('\t|\t')
            if name_class.replace('\t|\n','')=="scientific name":
                names[taxid] = name
    return names

def ibf_size_mb(args, bp, bins):
    return (math.ceil(-(1/((1-args.max_fp**(1/float(args.hash_functions)))**(1/float(args.hash_functions*(bp-args.kmer_size+1)))-1)))*bins)/8388608

def optimal_bins(nbins):
    return (math.floor(nbins/64)+1)*64

def get_rank_node(nodes, ranks, taxid, rank):
    t = taxid
    try:
        while ranks[t]!=rank and t!="1": t = nodes[t]
    except:
        return taxid
    
    return t if t!="1" else taxid

def bins_group(groups_len, fragment_size, overlap_length):
    group_nbins = {}
    # ignore overlap_length if too close to the size of the fragment size (*2) to avoid extreme numbers
    if fragment_size<=overlap_length*2: overlap_length=0
    for group, group_len in groups_len.items():
        # approximate extension in size with overlap_length (should be done by each sequence for the group)
        g_len = group_len + (math.floor(group_len/fragment_size)*overlap_length)
        group_nbins[group] = math.ceil(g_len/(fragment_size-overlap_length))
    return group_nbins

def estimate_bin_len(args, taxsbp_input_file, ncbi_nodes_file, use_assembly):
    ncbi_nodes, ncbi_ranks = read_nodes(ncbi_nodes_file)
    
    groups_len = defaultdict(int)
    for line in open(taxsbp_input_file, "r"):
        if use_assembly:
            _, seqlen, _, assembly = line.rstrip().split("\t")
            group = assembly
        else:
            fields = line.rstrip().split("\t") # slip by fields in case file has more
            seqlen = fields[1]
            taxid = fields[2]
            group = taxid if args.rank=="taxid" else get_rank_node(ncbi_nodes, ncbi_ranks, taxid, args.rank) # convert taxid into rank taxid

        sl = int(seqlen)
        if sl<args.kmer_size: continue
        groups_len[group]+= sl

    # number of groups and aux. variables
    ngroups = len(groups_len)
    min_group_len = min(groups_len.values())
    max_group_len = max(groups_len.values())
    sum_group_len = sum(groups_len.values())

    # special case with one group
    if ngroups==1:
        return math.floor(sum_group_len/64)

    # minimum number of bins possible (= number of groups) will generate a big and sparse IBF
    min_bins_optimal = optimal_bins(ngroups)
    # maximum number of bins possible (bin_length = min_group_len) will generate the smallest possible IBF
    max_bins_optimal = optimal_bins(sum(bins_group(groups_len, min_group_len, args.overlap_length).values()))
    # Min. possible size based on the maxium number of bins
    min_size_possible = ibf_size_mb(args, min_group_len, max_bins_optimal)

    if args.verbose:
        print_log("Group sizes for %s in bp: min %d, max %d, avg %d, sum %d\n" % (args.rank, min_group_len, max_group_len, math.ceil(sum_group_len/min_bins_optimal), sum_group_len))
        print_log("Minimum optimal number of bins %d (from %d bins)\n" % (min_bins_optimal,ngroups))
        print_log("Minimum estimated size %.2fMB in 1 bin\n" % ibf_size_mb(args, sum_group_len, 1))
        print_log("%d bins necessary to achieve %.2fMB (%dbp per bin)\n" % (max_bins_optimal, min_size_possible, min_group_len))
        print_log("%.2fMB necessary to achieve %d bins (%dbp per bin)\n" % (ibf_size_mb(args, max_group_len, min_bins_optimal), min_bins_optimal, max_group_len))

    # define maximum size based on an user input
    if args.max_bloom_size:
        if args.max_bloom_size < math.ceil(min_size_possible):
            print("--max-bloom-size is smaller than minimum possible %.2f. Try increasing --max-fp/--hash-functions" % min_size_possible)
            return 0
        max_bloom_size = args.max_bloom_size # user provided size limit
        min_nbins = min_bins_optimal # use the optimal
    else:
        max_bloom_size = 0 # no size limit
        min_nbins = min_bins_optimal*2 # define a threshold for stoping

    # Try to estimate best fit for size and number of bins
    bin_length = min_group_len  # start length on the minimun possible, max bins (biggest filter)
    bin_increment = math.ceil(sum_group_len/min_bins_optimal) # increment filter based on the avg.
    attempts=30 # number of iterations
    att_cont=0
    decrease_iter = 1 
    if args.verbose: print_log("\nmax_bloom_size %d, min_nbins %d\n" % (max_bloom_size, min_nbins))
    while att_cont<=attempts:
        att_cont+=1

        nbins = optimal_bins(sum(bins_group(groups_len, bin_length, args.overlap_length).values()))
        size_mb = ibf_size_mb(args, bin_length, nbins)
  
        if args.verbose: print_log("bin_length %d will generate %d bins with %.2fMB\n" % (bin_length, nbins, size_mb))

        if max_bloom_size: # has user input size limit
            if math.ceil(size_mb)>=max_bloom_size*0.99 and math.floor(size_mb)<=max_bloom_size*1.01: # if achieved in the goal range
                break
            elif size_mb>max_bloom_size: # index too big, step back and reduce increment
                decrease_iter=2 # on first time passing max_bloom_size start to decrease increment on every iteration
                bin_increment = -abs(math.floor(bin_increment/decrease_iter))
            elif size_mb<max_bloom_size:
                bin_increment = abs(math.floor(bin_increment/decrease_iter))
        else:
            if nbins>=min_nbins-64 and nbins<=min_nbins+64: # achieved low number of bins
                break
            elif nbins>min_nbins:
                bin_increment = abs(math.floor(bin_increment/decrease_iter))
            elif nbins<min_nbins:
                decrease_iter=2  # on first time below min_nbins start to decrease increment on every iteration
                bin_increment = -abs(math.floor(bin_increment/decrease_iter))

        if bin_increment==0: break

        bin_length+=bin_increment

        # if bin_length is growing bigger than the max group len. or smaller than overlap
        if bin_length>=max_group_len: 
            bin_length=max_group_len
            break
        elif bin_length<=args.overlap_length:
            bin_length=args.overlap_length+1
            break
    
    if args.verbose: print_log("-- %d bins -> %.2fMB (%dbp per bin) -- \n" % (nbins, size_mb , bin_length))

    return bin_length

def taxid_rank_up_to(taxid, filtered_nodes, fixed_ranks):
    if taxid!="0":
        original_rank = filtered_nodes[taxid][2]
        original_taxid = taxid
        while taxid!="0":
            if(filtered_nodes[taxid][2] in fixed_ranks):
                #everything below species (not being assembly) is counted as species+
                if "species+" in fixed_ranks and original_rank!="species" and original_rank!="assembly" and filtered_nodes[taxid][2]=="species":
                    return original_taxid, "species+"
                else:
                    return taxid, filtered_nodes[taxid][2]
            taxid = filtered_nodes[taxid][0]
        return "1", "root" #no standard rank identified
    else:
        return "0", ""

def print_final_report(reports, merged_filtered_nodes, seq_unc, total_reads, final_report_file, ranks):

    if not ranks:  
        all_ranks = False
        fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+','assembly']
    elif ranks[0]=="all":
        all_ranks = True  
        fixed_ranks = []
    else:
        all_ranks = False
        fixed_ranks = ['root'] + ranks

    # merge nodes
    all_filtered_nodes = {}
    for m in merged_filtered_nodes.values():
        all_filtered_nodes.update(m)

    # sum counts of each report
    merged_rep = defaultdict(int)
    for rep in reports.values():
        for leaf in rep.keys():
            merged_rep[leaf]+=rep[leaf]['count']

    final_rep = defaultdict(lambda: {'count': 0, 'rank': ""})
    
    # make cummulative sum of the counts on the lineage
    for leaf in merged_rep.keys():
        count = merged_rep[leaf]
        if all_ranks: # use all nodes on the tree
            t = leaf
            r = all_filtered_nodes[t][2]
        else: # use only nodes of the fixed ranks
            t, r = taxid_rank_up_to(leaf, all_filtered_nodes, fixed_ranks)

        while t!="0":
            final_rep[t]['count']+=count
            final_rep[t]['rank']=r
            if all_ranks:
                t = all_filtered_nodes[t][0]
                r = all_filtered_nodes[t][2] if t!="0" else ""
            else:
                t, r = taxid_rank_up_to(all_filtered_nodes[t][0], all_filtered_nodes, fixed_ranks)

    # build lineage after all entries were defined
    lineage = {}
    for assignment in final_rep.keys():
        lineage[assignment]=[]
        if all_ranks:
            t=assignment
            while t!="0":
                lineage[assignment].insert(0,t)
                t = all_filtered_nodes[t][0]
        else:
            t, r = taxid_rank_up_to(assignment, all_filtered_nodes, fixed_ranks)
            max_rank_idx = fixed_ranks.index(r) # get index of current rank
            while t!="0":
                # Add empty || if fixed rank is missing
                for i in range(max_rank_idx-fixed_ranks.index(r)):
                    lineage[assignment].insert(0,"")
                    max_rank_idx-=1

                lineage[assignment].insert(0,t)
                max_rank_idx-=1
                t, r = taxid_rank_up_to(all_filtered_nodes[t][0], all_filtered_nodes, fixed_ranks)

    with open(final_report_file, 'w') as frfile:
        frfile.write("unclassified" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ str(seq_unc) +"\t"+ str("%.5f" % ((seq_unc/total_reads)*100)) + "\n")
        if all_ranks:
            for assignment in sorted(lineage, key=lineage.get):
                frfile.write(final_rep[assignment]['rank'] +"\t"+ assignment +"\t"+ "|".join(lineage[assignment]) +"\t"+ all_filtered_nodes[assignment][1] +"\t"+ str(final_rep[assignment]['count']) +"\t"+ str("%.5f" % ((final_rep[assignment]['count']/total_reads)*100)) + "\n")
        else:
            for assignment in sorted(lineage, key=lambda k: (fixed_ranks.index(final_rep[k]['rank']), -final_rep[k]['count']), reverse=False):
                frfile.write(final_rep[assignment]['rank'] +"\t"+ assignment +"\t"+ "|".join(lineage[assignment]) +"\t"+ all_filtered_nodes[assignment][1] +"\t"+ str(final_rep[assignment]['count']) +"\t"+ str("%.5f" % ((final_rep[assignment]['count']/total_reads)*100)) + "\n")

def print_log(text):
    sys.stderr.write(text)
    sys.stderr.flush()

def get_output_line(file):
    readid=""
    classification=""
    kmercount=""
    done=False

    while True:
        try: # read first entry
            last_pos = file.tell() # save position in case of incomplete buffer
            line = file.readline() 
            if line[-1]!="\n": # incomplete line (buffer)
                file.seek(last_pos) # return the pointer to the beginning of the line
                raise ValueError # raise error

            readid, classification, kmercount = line.rstrip().split("\t") # may raise ValueError when it is the last line (added '\n')
            break
        except (ValueError, IndexError): # last line -> empty, no reads classified, broken line
            #print_log(file.name + " " + str(last_pos))
            #print_log(" breakline found (over)\n" if line=="\n" else " not over - " + line + "\n")
            # line=="\n" -> extra single  "\n" to tell the file is over (ganon-classify)
            if line=="\n": 
                done = True # do not start over
                break                
            else:
                time.sleep(1) # wait and try again
                continue

    return readid, classification, kmercount, done


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
        return 'Sequence({})'.format(', '.join(args))

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

    def write_bins(self, output_file):
        if self.bins:
            with open(output_file, 'w') as file:
                for line in self.bins:
                    file.write(line+"\n")

if __name__ == '__main__':
    main()
