#!/usr/bin/env python3
import argparse, sys, shutil
from ganon.util import *

class Config:

    version = '0.4.1'
    path_exec = {'build': "", 'classify': "", 'get_seq_info': ""}
    empty = False

    def __init__(self, which: str=None, **kwargs):

        parser = argparse.ArgumentParser(prog='ganon', description='ganon', conflict_handler='resolve')
        parser.add_argument('-v', '--version', action='version', version='version: %(prog)s ' + self.version, help="Show program's version number and exit.")

        build_parser = argparse.ArgumentParser(description='Build options', add_help=False)
        
        # Required
        build_group_required = build_parser.add_argument_group('required arguments')
        build_group_required.add_argument('-d', '--db-prefix',   type=str,            required=True,  help='Database output prefix (.ibf, .map, .tax, .gnn will be created)')
        build_group_required.add_argument('-i', '--input-files', type=str, nargs="*", required=False, help='Input reference sequence fasta files [.gz]')
        
        # Defaults
        build_group_optional = build_parser.add_argument_group('optional arguments')
        build_group_optional.add_argument('-r', '--rank',            type=str,            metavar='', default='species', help='Rank specific target for classification [species,genus,...]. use "leaves" to use the leaf taxonomic node assigned to each sequence as targets. If specified rank is not found in the lineage, use the leaf taxonomic node as target. Default: species')
        build_group_optional.add_argument('-s', '--specialization',  type=str,            metavar='', default="",        help='Add extra specialized "rank" as target for classification after taxonomic leaves. When selected --rank is set to leaves. Options: [sequence,file,assembly,custom]. "sequence" will use sequence accession as target. "file" uses the filename as target. "assembly" will use assembly info from NCBI as target. "custom" uses the 4th column of the file provided in --seq-info-file as target.')
        build_group_optional.add_argument('-k', '--kmer-size',       type=int,            metavar='', default=19,        help='The k-mer size for the interleaved bloom filter. Default: 19')
        build_group_optional.add_argument('-n', '--hash-functions',  type=int,            metavar='', default=3,         help='The number of hash functions for the interleaved bloom filter. Default: 3')
        build_group_optional.add_argument('-f', '--max-fp',          type=float,          metavar='', default=0.05,      help='Max. false positive rate for k-mer classification. Default: 0.05')
        build_group_optional.add_argument('-m', '--max-bloom-size',  type=int,            metavar='',                    help='Approx. maximum filter size in Megabytes (MB). Will estimate best --bin-length based on --kmer-size, --hash-functions and --max-fp  [Mutually exclusive --fixed-bloom-size]')
        build_group_optional.add_argument('-l', '--bin-length',      type=int,            metavar='',                    help='Maximum length (in bp) for each bin. Default: auto')
        build_group_optional.add_argument('-t', '--threads',         type=int,            metavar='', default=2,         help='Number of sub-processes/threads to use. Default: 2')
        build_group_optional.add_argument('--fixed-bloom-size',      type=int,            metavar='',                    help='Fixed size for filter in Megabytes (MB), will ignore --max-fp [Mutually exclusive --max-bloom-size] ')
        build_group_optional.add_argument('--fragment-length',       type=int,            metavar='', default=-1,        help='Fragment length (in bp). Set to 0 to not fragment sequences. Default: --bin-length - --overlap-length')
        build_group_optional.add_argument('--overlap-length',        type=int,            metavar='', default=300,       help='Fragment overlap length (in bp). Should be bigger than the read length used for classification. Default: 300')
        build_group_optional.add_argument('--seq-info-mode',         type=str, nargs="*", metavar='', default=["auto"],  help='Automatic mode to retrieve tax. info and seq. length. [auto,eutils] or one or more accession2taxid files from NCBI [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot]. auto will either use eutils for less than 50000 input sequences or nucl_gb nucl_wgs. Alternatively a file can be directly provided (see --seq-info-file). Default: auto')
        build_group_optional.add_argument('--seq-info-file',         type=str,            metavar='',                    help='Tab-separated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> specialization]) [Mutually exclusive --seq-info-mode]')
        build_group_optional.add_argument('--taxdump-file',          type=str, nargs="*", metavar='',                    help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
        build_group_optional.add_argument('--input-directory',       type=str,            metavar='', default="",        help='Directory containing input files')
        build_group_optional.add_argument('--input-extension',       type=str,            metavar='', default="",        help='Extension of files to use with --input-directory (provide it without * expansion, e.g. ".fna.gz")')
        build_group_optional.add_argument('--write-seq-info-file',   action='store_true',                                help='Write sequence information to DB_PREFIX.seqinfo.txt')
        build_group_optional.add_argument('--verbose',               action='store_true',                                help='Verbose output mode')
        build_group_optional.add_argument('--quiet',                 action='store_true',                                help='Quiet output mode')
        build_group_optional.add_argument('--ganon-path',            type=str,            metavar='', default="",        help=argparse.SUPPRESS)
        build_group_optional.add_argument('--n-refs',                type=int,            metavar='',                    help=argparse.SUPPRESS)
        build_group_optional.add_argument('--n-batches',             type=int,            metavar='',                    help=argparse.SUPPRESS)
        
        ####################################################################################################

        update_parser = argparse.ArgumentParser(description='Update options', add_help=False)

        # Required
        update_group_required = update_parser.add_argument_group('required arguments')
        update_group_required.add_argument('-d', '--db-prefix',   type=str,            required=True,  help='Database input prefix (.ibf, .map, .tax, .gnn)')
        update_group_required.add_argument('-i', '--input-files', type=str, nargs="*", required=False, help='Input reference sequence fasta files [.gz] to be included to the database. Complete set of updated sequences should be provided when using --update-complete')
        
        # Defaults
        update_group_optional = update_parser.add_argument_group('optional arguments')
        update_group_optional.add_argument('-o', '--output-db-prefix', type=str,            metavar='',                   help='Output database prefix (.ibf, .map, .tax, .gnn). Default: overwrite current --db-prefix')
        update_group_optional.add_argument('-t', '--threads',          type=int,            metavar='', default=2,        help='Number of sub-processes/threads to use. Default: 2')
        update_group_optional.add_argument('-s', '--specialization',   type=str,            metavar='', default="",       help='Change specialization mode. Can only be used if database was built with some specialization. Options: [sequence,file,assembly,custom]. "sequence" will use sequence accession as target. "file" uses the filename as target. "assembly" will use assembly info from NCBI as target. "custom" uses the 4th column of the file provided in --seq-info-file as target.')
        update_group_optional.add_argument('--seq-info-mode',          type=str, nargs="*", metavar='', default=["auto"], help='Automatic mode to retrieve tax. info and seq. length. [auto,eutils] or one or more accession2taxid files from NCBI [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot]. auto will either use eutils for less than 50000 input sequences or nucl_gb nucl_wgs. Alternatively a file can be directly provided (see --seq-info-file). Default: auto')
        update_group_optional.add_argument('--seq-info-file',          type=str,            metavar='',                   help='Tab-separated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly id]) [Mutually exclusive --seq-info]')
        update_group_optional.add_argument('--taxdump-file',           type=str, nargs="*", metavar='',                   help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
        update_group_optional.add_argument('--input-directory',        type=str,            metavar='', default="",       help='Directory containing input files')
        update_group_optional.add_argument('--input-extension',        type=str,            metavar='', default="",       help='Extension of files to use with --input-directory (provide it without * expansion, e.g. ".fna.gz")')
        update_group_optional.add_argument('--update-complete',        action='store_true',                               help='Update adding and removing sequences. Input files should represent the complete updated set of references, not only new sequences.')
        update_group_optional.add_argument('--write-seq-info-file',    action='store_true',                               help='Write sequence information to DB_PREFIX.seqinfo.txt')
        update_group_optional.add_argument('--verbose',                action='store_true',                               help='Verbose output mode')
        update_group_optional.add_argument('--quiet',                  action='store_true',                               help='Quiet output mode')
        update_group_optional.add_argument('--ganon-path',             type=str,            metavar='', default="",       help=argparse.SUPPRESS)
        update_group_optional.add_argument('--n-refs',                 type=int,            metavar='',                   help=argparse.SUPPRESS)
        update_group_optional.add_argument('--n-batches',              type=int,            metavar='',                   help=argparse.SUPPRESS)

        ####################################################################################################

        classify_parser = argparse.ArgumentParser(description='Classification options', add_help=False)

        # Required
        classify_group_required = classify_parser.add_argument_group('required arguments')
        classify_group_required.add_argument('-d', '--db-prefix',    type=str, nargs="*", required=True,                                             help='Database input prefix[es]')
        classify_group_required.add_argument('-s', '--single-reads', type=str, nargs="*", required=False, metavar='reads.fq[.gz]',                   help='Multi-fastq[.gz] file[s] to classify')
        classify_group_required.add_argument('-p', '--paired-reads', type=str, nargs="*", required=False, metavar='reads.1.fq[.gz] reads.2.fq[.gz]', help='Multi-fastq[.gz] pairs of file[s] to classify')

        # Defaults
        classify_group_optional = classify_parser.add_argument_group('optional arguments')
        classify_group_optional.add_argument('-o', '--output-prefix',       type=str,              metavar='', help='Output prefix for .lca and .rep. Empty to output to STDOUT (only .lca will be printed)')
        classify_group_optional.add_argument('-c', '--hierarchy-labels',    type=str,   nargs="*", metavar='', help='Hierarchy definition, one for each database input. Can also be a string, but input will be sorted to define order (e.g. 1 1 2 3). Default: 1')
        classify_group_optional.add_argument('-k', '--min-kmers',           type=float, nargs="*", metavar='', help='Min. percentage of k-mers matching to consider a read assigned. Single value or one per database (e.g. 0.5 0.7 1 0.25). Default: 0.25 [Mutually exclusive --max-error]')
        classify_group_optional.add_argument('-e', '--max-error',           type=int,   nargs="*", metavar='', help='Max. number of errors allowed. Single value or one per database (e.g. 3 3 4 0) [Mutually exclusive --min-kmers]')
        classify_group_optional.add_argument('-u', '--max-error-unique',    type=int,   nargs="*", metavar='', help='Max. number of errors allowed for unique assignments after filtering. Matches below this error rate will not be discarded, but assigned to a parent taxonomic level. Single value or one per hierarchy (e.g. 0 1 2). -1 to disable. Default: -1')    
        classify_group_optional.add_argument('-l', '--strata-filter',       type=int,   nargs="*", metavar='', help='Additional errors allowed (relative to the best match) to filter and select matches. Single value or one per hierarchy (e.g. 0 1 2). -1 to disable filtering. Default: 0')    
        classify_group_optional.add_argument('-f', '--offset',              type=int,              metavar='', help='Number of k-mers to skip during classification. Can speed up analysis but may reduce recall. (e.g. 1 = all k-mers, 3 = every 3rd k-mer). Default: 2')    
        classify_group_optional.add_argument('-t', '--threads',             type=int,              metavar='', help='Number of sub-processes/threads to use. Default: 3')
        classify_group_optional.add_argument('-r', '--ranks',               type=str,   nargs="*", metavar='', help='Ranks to show in the report (.tre). "all" for all identified ranks. empty for default ranks: superkingdom phylum class order family genus species assembly. This file can be re-generated with the ganon report command.')
        classify_group_optional.add_argument('--output-all',                action='store_true',               help='Output an additional file with all matches (.all). File can be very large.')
        classify_group_optional.add_argument('--output-unclassified',       action='store_true',               help='Output an additional file with unclassified read headers (.unc)')
        classify_group_optional.add_argument('--output-single',             action='store_true',               help='When using multiple hierarchical levels, output everything in one file instead of one per hierarchy')
        classify_group_optional.add_argument('--verbose',                   action='store_true',               help='Verbose output mode')
        classify_group_optional.add_argument('--quiet',                     action='store_true',               help='Quiet output mode')
        classify_group_optional.add_argument('--ganon-path',                type=str, default="",  metavar='', help=argparse.SUPPRESS) 
        classify_group_optional.add_argument('--n-reads',                   type=int,              metavar='', help=argparse.SUPPRESS)
        classify_group_optional.add_argument('--n-batches',                 type=int,              metavar='', help=argparse.SUPPRESS)
        
        ####################################################################################################
        
        report_parser = argparse.ArgumentParser(description='Report options', add_help=False)

        # Required
        report_group_required = report_parser.add_argument_group('required arguments')
        report_group_required.add_argument('-i', '--rep-files',     type=str, nargs="*", required=False, help='One or more *.rep files from ganon classify')
        report_group_required.add_argument('-o', '--output-prefix', type=str,            required=True,  help='Output prefix for report file "{output_prefix}.tre". In case of multiple files, the base input filename will be appended at the end of the output file "{output_prefix + FILENAME}.tre"')
    
        # Defaults
        report_group_optional = report_parser.add_argument_group('optional arguments')
        report_group_optional.add_argument('-d', '--db-prefix',      type=str, nargs="*", metavar='', default=[],      help='Database prefix[es] used for classification (in any order). Only ".tax" file is required. If not provided, new taxonomy will be downloaded')
        report_group_optional.add_argument('-f', '--output-format',  type=str,            metavar='', default="tsv",   help='Output format [text, tsv, csv]. text outputs a tabulated formatted text file for better visualization. Default: tsv')
        report_group_optional.add_argument('-e', '--report-type',    type=str,            metavar='', default="reads", help='Type of report to generate [reads, matches]. Default: reads')
        report_group_optional.add_argument('-r', '--ranks',          type=str, nargs="*", metavar='', default=[],      help='Ranks to report ["", "all", custom list] "all" for all possible ranks. empty for default ranks (superkingdom phylum class order family genus species assembly). Default: ""')
        report_group_optional.add_argument('-s', '--sort',           type=str,            metavar='', default="",      help='Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage (with --ranks all)')
        report_group_optional.add_argument('-y', '--split-hierarchy',action='store_true',                              help='Split output reports by hierarchy (from ganon classify --hierarchy-labels). If activated, the output files will be named as "{output_prefix}.{hierarchy}.tre"')
        report_group_optional.add_argument('-p', '--skip-hierarchy', type=str, nargs="*", metavar='', default=[],      help='One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels)')
        report_group_optional.add_argument('-k', '--keep-hierarchy', type=str, nargs="*", metavar='', default=[],      help='One or more hierarchies to keep in the report (from ganon classify --hierarchy-labels)')
        report_group_optional.add_argument('--taxdump-file',         type=str, nargs="*", metavar='', default=[],      help='Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)')
        report_group_optional.add_argument('--input-directory',      type=str,            metavar='', default="",      help='Directory containing input files')
        report_group_optional.add_argument('--input-extension',      type=str,            metavar='', default="",      help='Extension of files to use with --input-directory (provide it without * expansion, e.g. ".rep")')
        report_group_optional.add_argument('--verbose',              action='store_true',             default=False,   help='Verbose output mode')
        report_group_optional.add_argument('--quiet',                action='store_true',             default=False,   help='Quiet output mode')

        ####################################################################################################

        table_parser = argparse.ArgumentParser(description='Table options', add_help=False)

        # Required
        table_group_required = table_parser.add_argument_group('required arguments')
        table_group_required.add_argument('-i', '--tre-files',   type=str, nargs="*", required=False, help='Report files (.tre) from ganon classify/report to make the table')
        table_group_required.add_argument('-o', '--output-file', type=str,            required=True,  help='Output filename for the table')

        # Defaults
        table_group_optional = table_parser.add_argument_group('optional arguments')
        table_group_optional.add_argument('-l', '--output-value',             type=str,   metavar='', default="counts",     help="Output value on the table [percentage, counts]. percentage values are reported between [0-1]. Default: counts")
        table_group_optional.add_argument('-f', '--output-format',            type=str,   metavar='', default="tsv",        help='Output format [tsv, csv]. Default: tsv')
        table_group_optional.add_argument('-t', '--top-sample',               type=int,   metavar='', default=0,            help="Top hits of each sample individually")
        table_group_optional.add_argument('-a', '--top-all',                  type=int,   metavar='', default=0,            help="Top hits of all samples (ranked by percentage)")
        table_group_optional.add_argument('-r', '--rank',                     type=str,   metavar='', default=None,         help="Define specific rank to report. Empty will report all ranks (only direct read assignments - col. 6 from the report files)")
        table_group_optional.add_argument('-m', '--min-occurrence',           type=int,   metavar='', default=0,            help="# occurrence of a taxa among reports to be kept [1-*]")
        table_group_optional.add_argument('-p', '--min-occurrence-percentage',type=float, metavar='', default=0,            help="%% occurrence of a taxa among reports to be kept [0-1]")
        table_group_optional.add_argument('--header',                         type=str,   metavar='', default="name",       help='Header information [name, taxid, lineage]. Default: name')
        table_group_optional.add_argument('--add-unclassified',               action='store_true',    default=False,        help="Add column with unclassified count/percentage")
        table_group_optional.add_argument('--add-filtered',                   action='store_true',    default=False,        help="Add column with filtered count/percentage")
        table_group_optional.add_argument('--skip-zeros',                     action='store_true',    default=False,        help="Do not print lines with only zero count/percentage")
        table_group_optional.add_argument('--transpose',                      action='store_true',    default=False,        help="Transpose output table (taxa as cols and files as rows)")
        table_group_optional.add_argument('--input-directory',                type=str,  metavar='',  default="",           help='Directory containing input files')
        table_group_optional.add_argument('--input-extension',                type=str,  metavar='',  default="",           help='Extension of files to use with --input-directory (provide it without * expansion, e.g. ".tre")')
        table_group_optional.add_argument('--verbose',                        action='store_true',    default=False,        help='Verbose output mode')
        table_group_optional.add_argument('--quiet',                          action='store_true',    default=False,        help='Quiet output mode')

        ####################################################################################################

        filter_parser = argparse.ArgumentParser(description='Table options', add_help=False)
        filter_arguments = filter_parser.add_argument_group('filter arguments')
        filter_arguments.add_argument('--min-count',      type=int,            metavar='', default=0,  help="Minimum number of counts to keep the taxa")
        filter_arguments.add_argument('--min-percentage', type=float,          metavar='', default=0,  help="Minimum percentage of counts to keep the taxa [0-1]")
        filter_arguments.add_argument('--names',          type=str, nargs="*", metavar='', default=[], help="Show only entries matching exact names of the provided list")
        filter_arguments.add_argument('--names-with',     type=str, nargs="*", metavar='', default=[], help="Show entries containing full or partial names of the provided list")
        filter_arguments.add_argument('--taxids',         type=str, nargs="*", metavar='', default=[], help='One or more taxids to report (including children taxa)')
        
        subparsers = parser.add_subparsers()
        
        build = subparsers.add_parser('build', help='Build ganon database', parents=[build_parser])
        build.set_defaults(which='build')

        update = subparsers.add_parser('update', help='Update ganon database', parents=[update_parser])
        update.set_defaults(which='update')

        classify = subparsers.add_parser('classify', help='Classify reads', parents=[classify_parser])
        classify.set_defaults(which='classify')

        report = subparsers.add_parser('report', help='Generate reports', parents=[filter_parser,report_parser])
        report.set_defaults(which='report')

        table = subparsers.add_parser('table', help='Generate table from reports', parents=[filter_parser,table_parser])
        table.set_defaults(which='table')

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

        seq_info_mode_options = ["auto","eutils","nucl_gb","nucl_wgs","nucl_est","nucl_gss","pdb","prot","dead_nucl","dead_wgs","dead_prot"]

        if self.empty is True:
            print_log("Please provide one or more arguments")
            return False

        if self.verbose is True:
            self.quiet = False
        elif self.quiet is True:
            self.verbose = False

        if self.which == 'build':
            if not check_taxdump(self.taxdump_file):
                return False

            if not check_input_directory(self.input_files, self.input_directory, self.input_extension):
                return False

            if self.fragment_length > 0 and self.overlap_length > self.fragment_length:
                print_log("--overlap-length cannot be bigger than --fragment-length")
                return False

            if self.fixed_bloom_size and not self.bin_length:
                print_log("please set the --bin-length to use --fixed-bloom-size")
                return False

            if self.max_fp <= 0:
                print_log("--max-fp has to be bigger than 0")
                return False

            if self.specialization:
                valid_spec = self.validate_specialization()
                if not valid_spec:
                    return False

            if not all([sim in seq_info_mode_options for sim in self.seq_info_mode]):
                print_log("Invalid --seq-info-mode. Options: " + " ".join(seq_info_mode_options))
                return False

        elif self.which=='update':
            if not check_taxdump(self.taxdump_file):
                return False

            if not check_input_directory(self.input_files, self.input_directory, self.input_extension):
                return False

            if not check_db(self.db_prefix):
                return False

            if self.specialization: 
                valid_spec = self.validate_specialization()
                if not valid_spec: return False

            if not all([sim in seq_info_mode_options for sim in self.seq_info_mode]):
                print_log("Invalid --seq-info-mode. Options: " + " ".join(seq_info_mode_options))
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

            if self.skip_hierarchy and self.keep_hierarchy:
                print_log("--skip-hierarchy and --keep-hierarchy are mutually exclusive")
                return False

            if not check_taxdump(self.taxdump_file):
                return False

            if not check_input_directory(self.rep_files, self.input_directory, self.input_extension):
                return False

            if self.db_prefix:
                dbp=[]
                # add ".tax" if was passed as prefix
                for prefix in self.db_prefix:
                    if prefix.endswith(".tax"):
                        dbp.append(prefix)
                    else:
                        dbp.append(prefix+".tax")
                # check if files exists
                dbp_ok = check_files(dbp)
                # report files not found and stop
                if len(dbp_ok)!=len(dbp):
                    print_log(",".join(set(dbp).difference(dbp_ok)))
                    return False

        elif self.which=='table':
            if not check_input_directory(self.tre_files, self.input_directory, self.input_extension):
                return False

            if self.min_occurrence < 0:
                print_log("Invalid value for --min-occurrence (>0)")
                return False

            if self.min_occurrence_percentage < 0 or self.min_occurrence_percentage > 1:
                print_log("Invalid value for --min-occurrence-percentage [0-1]")
                return False

        return True

    def validate_specialization(self):
        if self.specialization not in ["sequence","file","assembly","custom"]:
            print_log("Invalid value for --specialization")
            return False
        
        if self.specialization=="custom" and not self.seq_info_file:
            print_log("--seq-info-file should be provided with --specialization custom")
            return False

        if self.specialization!="custom" and self.seq_info_file:
            print_log("--specialization custom is required with --seq-info-file")
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

            ganon_get_seq_info_paths = [self.ganon_path, self.ganon_path+"scripts/", self.ganon_path+"../scripts/"] if self.ganon_path else [None, "scripts/"]
            for p in ganon_get_seq_info_paths:
                self.path_exec['get_seq_info'] = shutil.which("ganon-get-seq-info.sh", path=p)
                if self.path_exec['get_seq_info'] is not None: break
            if self.path_exec['get_seq_info'] is None:
                print_log("ganon-get-seq-info.sh script was not found. Please inform a specific path with --ganon-path")
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

def check_input_directory(input_files, input_directory, input_extension):         
    if not input_files and not input_directory:
        print_log("Please provide files with --input-files and/or --input-directory with --input-extension")
        return False
    elif input_directory and not input_extension:
        print_log("Please provide the --input-extension when using --input-directory")
        return False
    elif input_directory and "*" in input_extension:
        print_log("Please do not use wildcards (*) in the --input-extension")
        return False
    return True

def check_taxdump(taxdump_file):
    if taxdump_file:
        if ((len(taxdump_file)==1 and not taxdump_file[0].endswith(".tar.gz")) or len(taxdump_file)>3):
            print_log("Please provide --taxdump-file taxdump.tar.gz OR --taxdump-file nodes.dmp names.dmp OR --taxdump-file nodes.dmp names.dmp merged.dmp OR leave it empty for automatic download")
            return False
        else:
            for f in taxdump_file:
                if not os.path.isfile(f):
                    print_log("File not found: " + f)
                    return False
    return True

def check_db(prefix):
    for db_file_type in [".ibf", ".map", ".tax", ".gnn"]:
        if not os.path.isfile(prefix+db_file_type):
            print_log("Incomplete database [" + prefix  + "] (.ibf, .map, .tax and .gnn)")
            return False
    return True
