#!/usr/bin/env python3
import argparse, sys, shutil
from ganon.util import *

class Config:

    version = '0.3.3'
    path_exec = {'build': "", 'classify': "", 'get_len_taxid': ""}
    empty = False

    def __init__(self, which: str=None, **kwargs):

        parser = argparse.ArgumentParser(prog='ganon', description='ganon', conflict_handler='resolve')
        parser.add_argument('-v', '--version', action='version', version='version: %(prog)s ' + self.version, help="Show program's version number and exit.")

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
        build_group_optional.add_argument('--write-seq-info-file', default=False, action='store_true', help='Write sequence information to DB_PREFIX.seqinfo.txt')
        build_group_optional.add_argument('--verbose', default=False, action='store_true', help='Verbose output mode')
        build_group_optional.add_argument('--quiet', default=False, action='store_true', help='Quiet output mode (only errors and warnings to the stderr)')
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
        update_group_optional.add_argument('--write-seq-info-file', default=False, action='store_true', help='Write sequence information to DB_PREFIX.seqinfo.txt')
        update_group_optional.add_argument('--verbose', default=False, action='store_true', help='Verbose output mode')
        update_group_optional.add_argument('--quiet', default=False, action='store_true', help='Quiet output mode (only errors and warnings to the stderr)')
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
        classify_group_optional.add_argument('--verbose', default=False, action='store_true',  help='Verbose output mode')
        classify_group_optional.add_argument('--quiet', default=False, action='store_true', help='Quiet output mode (only errors and warnings to the stderr)')
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
        report_group_optional.add_argument('--verbose', default=False, action='store_true',  help='Verbose output mode')
        report_group_optional.add_argument('--quiet', default=False, action='store_true', help='Quiet output mode (only errors and warnings to the stderr)')
        
        ####################################################################################################

        subparsers = parser.add_subparsers()
        
        build = subparsers.add_parser('build', help='Build ganon database', parents=[build_parser])
        build.set_defaults(which='build')

        update = subparsers.add_parser('update', help='Update ganon database', parents=[update_parser])
        update.set_defaults(which='update')

        classify = subparsers.add_parser('classify', help='Classify reads', parents=[classify_parser])
        classify.set_defaults(which='classify')

        report = subparsers.add_parser('report', help='Generate reports', parents=[report_parser])
        report.set_defaults(which='report')

        # Passing arguments internally from call main(which, **kwargs)
        if which is not None:
            # Set which as the first parameter (mandatory)
            list_kwargs = [which]
            for arg,value in kwargs.items():
                # convert all others to argparse format (eg: input_files to --input-files)
                arg_formatted = "--" + arg.replace('_', '-')
                if isinstance(value, list): # unpack if list
                    list_kwargs.append(arg_formatted)
                    list_kwargs.extend(value)
                elif type(value)==bool and value is True: # add only arg if boolean flag activated
                    list_kwargs.append(arg_formatted)
                elif value:
                    list_kwargs.append(arg_formatted)
                    list_kwargs.append(str(value))
            # Parse from list saving arguments to this class
            parser.parse_args(list_kwargs, namespace=self)
        else:
            # parse from default CLI sys.argv saving arguments to this class
            parser.parse_args(namespace=self) 
            if len(sys.argv)==1: 
                parser.print_help()
                self.empty = True

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'Config({})'.format(', '.join(args))

    def validate(self):
        if self.empty is True:
            print_log("Please provide one or more arguments")
            return False

        if self.verbose is True: 
            self.quiet=False
        elif self.quiet is True:
            self.verbose=False
            
        if self.which in ['build','update']:
            if self.taxdump_file and ((len(self.taxdump_file)==1 and not self.taxdump_file[0].endswith(".tar.gz")) or len(self.taxdump_file)>3):
                print_log("Please provide --taxdump-file taxdump.tar.gz or --taxdump-file nodes.dmp names.dmp [merged.dmp] or leave it empty for automatic download")
                return False

            if not self.input_files and not self.input_directory:
                print_log("Please provide files with --input-files and/or --input-directory with --input-extension")
                return False
            elif self.input_directory and not self.input_extension:
                print_log("Please provide the --input-extension when using --input-directory")
                return False
            elif self.input_directory and "*" in self.input_extension:
                print_log("Please do not use wildcards (*) in the --input-extension")
                return False

            if self.which=='update':
                if not check_db(self.db_prefix):
                    return False

            if self.which=='build': #If set (!=0), should be smaller than fragment
                if self.fragment_length>0 and self.overlap_length > self.fragment_length:
                    print_log("--overlap-length cannot be bigger than --fragment-length")
                    return False

                if self.fixed_bloom_size and not self.bin_length:
                    print_log("please set the --bin-length to use --fixed-bloom-size")
                    return False

                if self.max_fp<=0:
                    print_log("--max-fp has to be bigger than 0")
                    return False
            
        elif self.which=='classify':
            for prefix in self.db_prefix:
                if not check_db(prefix):
                    return False

            if not self.single_reads and not self.paired_reads:
                print_log("Please provide file[s] with --single-reads or --paired-reads")
                return False

            len_single_reads = 0
            if self.single_reads: 
                self.single_reads = check_files(self.single_reads)
                len_single_reads = len(self.single_reads)
            len_paired_reads = 0
            if self.paired_reads: 
                self.paired_reads = check_files(self.paired_reads)
                len_paired_reads = len(self.paired_reads)
            
            if len_paired_reads % 2 != 0:
                print_log("Invalid paired reads")
                return False

            if len_single_reads+len_paired_reads==0:
                print_log("No valid input files to classify")
                return False

        elif self.which=='report':
            for prefix in self.db_prefix:
                if not check_db(prefix):
                    return False

            if not os.path.isfile(self.rep_file):
                print_log("File not found [" + self.rep_file + "]")
                return False

        return True

    def set_paths(self):
        missing_path = False
        if self.which in ['build','update']:
            self.ganon_path = self.ganon_path + "/" if self.ganon_path else ""

            # if path is given, look for binaries only there
            ganon_build_paths = [self.ganon_path, self.ganon_path+"build/"] if self.ganon_path else [None, "build/"]
            for p in ganon_build_paths:
                self.path_exec['build'] = shutil.which("ganon-build", path=p)
                if self.path_exec['build'] is not None: break
            if self.path_exec['build'] is None:
                print_log("ganon-build binary was not found. Please inform a specific path with --ganon-path")
                missing_path = True

            ganon_get_len_taxid_paths = [self.ganon_path, self.ganon_path+"scripts/", self.ganon_path+"../scripts/"] if self.ganon_path else [None, "scripts/"]
            for p in ganon_get_len_taxid_paths:
                self.path_exec['get_len_taxid'] = shutil.which("ganon-get-len-taxid.sh", path=p)
                if self.path_exec['get_len_taxid'] is not None: break
            if self.path_exec['get_len_taxid'] is None:
                print_log("ganon-get-len-taxid.sh script was not found. Please inform a specific path with --ganon-path")
                missing_path = True

        elif self.which in ['classify']:
            self.ganon_path = self.ganon_path + "/" if self.ganon_path else ""

            ganon_classify_paths = [self.ganon_path, self.ganon_path+"build/"] if self.ganon_path else [None, "build/"]
            for p in ganon_classify_paths:
                self.path_exec['classify'] = shutil.which("ganon-classify", path=p)
                if self.path_exec['classify'] is not None: break
            if self.path_exec['classify'] is None:
                print_log("ganon-classify binary was not found. Please inform a specific path with --ganon-path")
                missing_path = True

        return True if not missing_path else False