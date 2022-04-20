#!/usr/bin/env python3

import argparse
import sys
import shutil
from ganon.util import *


class Config:

    version = "1.2.0"
    path_exec = {"build": "", "classify": "", "get_seq_info": "", "genome_updater": ""}
    empty = False

    choices_taxonomy = ["ncbi", "gtdb", "none"] # get from multitax
    choices_og = ["archaea", "bacteria", "fungi", "human", "invertebrate", "metagenomes", "other", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"]
    choices_db_source = ["refseq", "genbank"]
    choices_level = ["assembly", "name", "custom"]
    choices_ncbi_sequence_info = ["eutils","nucl_gb","nucl_wgs","nucl_est","nucl_gss","pdb","prot","dead_nucl","dead_wgs","dead_prot"]
    choices_ncbi_file_info = ["refseq", "genbank"]

    def __init__(self, which: str=None, **kwargs):

        parser = argparse.ArgumentParser(prog="ganon", description="ganon", conflict_handler="resolve")
        parser.add_argument("-v", "--version", action="version", version="version: %(prog)s " + self.version, help="Show program's version number and exit.")

        ####################################################################################################

        build_default_parser = argparse.ArgumentParser(add_help=False)

        build_default_required_args = build_default_parser.add_argument_group("required arguments")
        build_default_required_args.add_argument("-d", "--db-prefix", type=str, required=True, help="Database output prefix")

        build_default_important_args = build_default_parser.add_argument_group("important arguments")
        build_default_important_args.add_argument("-x", "--taxonomy", type=str, metavar="", default="ncbi", help="Downloads and build taxonomy tree to enable taxonomic classification and reports [" + ",".join(self.choices_taxonomy) + "]", choices=self.choices_taxonomy)
        build_default_important_args.add_argument("-t", "--threads",  type=int, metavar="", default=2,      help="")

        build_default_advanced_args = build_default_parser.add_argument_group("advanced arguments")
        build_default_advanced_args.add_argument("-p", "--max-fp",         type=int_or_float(0,1), metavar="", default=0.05, help="Max. false positive rate for bloom filters [Mutually exclusive --filter-size].")
        build_default_advanced_args.add_argument("-f", "--filter-size",    type=float,             metavar="", default=0,    help="Fixed size for filter in Megabytes (MB) [Mutually exclusive --max-fp]")
        build_default_advanced_args.add_argument("-k", "--kmer-size",      type=int,               metavar="", default=19,   help="The k-mer size to split sequences.")
        build_default_advanced_args.add_argument("-w", "--window-size",    type=int,               metavar="", default=32,   help="The window-size to build filter with minimizers.")
        build_default_advanced_args.add_argument("-s", "--hash-functions", type=int,               metavar="", default=0,    help="The number of hash functions for the interleaved bloom filter [0-5]. 0 to detect optimal value.", choices=range(6))

        build_default_other_args = build_default_parser.add_argument_group("options")
        build_default_other_args.add_argument("--restart",    action="store_true",                         help="Restart build from scratch even if files are found. Will overwrite {db_prefix}.tax, {db_prefix}.ibf) and {db_prefix}_tmp/")
        build_default_other_args.add_argument("--verbose",    action="store_true",                         help="Verbose output mode")
        build_default_other_args.add_argument("--quiet",      action="store_true",                         help="Quiet output mode")
        build_default_other_args.add_argument("--ganon-path", type=str,            metavar="", default="", help=argparse.SUPPRESS)
        build_default_other_args.add_argument("--n-refs",     type=int,            metavar="",             help=argparse.SUPPRESS)
        build_default_other_args.add_argument("--n-batches",  type=int,            metavar="",             help=argparse.SUPPRESS)

        ####################################################################################################

        build_parser = argparse.ArgumentParser(add_help=False)

        build_required_args = build_parser.add_argument_group("required arguments")
        build_required_args.add_argument("-g", "--organism-group", type=str, nargs="*", metavar="", required=True, help="One or more organims groups [" + ",".join(self.choices_og) + "]", choices=self.choices_og)

        build_download_args = build_parser.add_argument_group("download arguments")
        build_download_args.add_argument("-b", "--source",           type=str,            nargs="*", default="refseq", metavar="", help="[" + ",".join(self.choices_db_source) + "]", choices=self.choices_db_source)
        build_download_args.add_argument("-c", "--complete-genomes", action="store_true",                                          help="Download only sub-set of complete genomes")
        build_download_args.add_argument("-o", "--top",              type=int,                       default=0,        metavar="", help="Download top organims for each taxa (only possible for --taxonomy ncbi). 0 for all.")
        build_download_args.add_argument("-u", "--genome-updater",   type=str,                                         metavar="", help="Additional genome_updater parameters (https://github.com/pirovc/genome_updater)")

        ####################################################################################################

        build_custom_parser = argparse.ArgumentParser(add_help=False)

        build_custom_required_args = build_custom_parser.add_argument_group("required arguments")
        build_custom_required_args.add_argument("-i", "--input",           type=str,          nargs="*",        metavar="", help="Input file(s) and/or folder(s) [Mutually exclusive --input-file]")
        build_custom_required_args.add_argument("-e", "--input-extension", type=str,          default="fna.gz", metavar="", help="Required if --input contains folder(s). Wildcards are not supported (e.g. *).")

        build_custom_args = build_custom_parser.add_argument_group("custom arguments")
        build_custom_args.add_argument("-n", "--input-file",           type=file_exists,                                  metavar="", help="Tab separated file with information for each file and target: target <tab> tax.node <tab> specialization <tab> file [Mutually exclusive --input]")
        build_custom_args.add_argument("-a", "--input-target",         type=str,                                          metavar="", help="Target to use [file, sequence]. By default: 'file' if multiple input files are provided, 'sequence' if a single file is provided. Using 'file' is recommended and will speed-up the building process", choices=["file", "sequence"])
        build_custom_args.add_argument("-l", "--level",                type=str,                                          metavar="", help="Use a specialized target to build the database. Options: any available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy). Further specialization options [" + ",".join(self.choices_level) + "]. assembly will retrieve and use the ncbi assembly accession. name will retrieve and use the organism name. custom requires and uses the specialization field in the --input-file. By default, last level is set to the --input-target")
        build_custom_args.add_argument("-m", "--taxonomy-files",       type=str, nargs="*",                               metavar="", help="Specific files for taxonomy - otherwise files will be downloaded")

        ncbi_args = build_custom_parser.add_argument_group("ncbi arguments")
        ncbi_args.add_argument("-r", "--ncbi-sequence-info",    type=str, nargs="*", default=[],                          metavar="", help="Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information. [" + ",".join(self.choices_ncbi_sequence_info) + " or one or more accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/]. By default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise.")
        ncbi_args.add_argument("-q", "--ncbi-file-info",        type=str, nargs="*", default=self.choices_ncbi_file_info, metavar="", help="Downloads assembly_summary files to extract target information. [" + ",".join(self.choices_ncbi_file_info) + " or one or more assembly_summary files from https://ftp.ncbi.nlm.nih.gov/genomes/]")          

        ####################################################################################################

        update_custom_parser = argparse.ArgumentParser(add_help=False)
        update_custom_parser.add_argument("-d", "--db-prefix",   type=str,            required=True,  help="Database input prefix")
        
        ####################################################################################################

        update_parser = argparse.ArgumentParser(add_help=False)

        # Required
        update_group_required = update_parser.add_argument_group("required arguments")
        update_group_required.add_argument("-d", "--db-prefix",   type=str,            required=True,  help="Database input prefix (.ibf, .map, .tax, .gnn)")
        update_group_required.add_argument("-i", "--input-files", type=str, nargs="*", required=False, help="Input reference sequence fasta files [.gz] to be included to the database. Complete set of updated sequences should be provided when using --update-complete")

        # Defaults
        update_group_optional = update_parser.add_argument_group("optional arguments")
        update_group_optional.add_argument("-o", "--output-db-prefix", type=str,            metavar="",                   help="Output database prefix (.ibf, .map, .tax, .gnn). Default: overwrite current --db-prefix")
        update_group_optional.add_argument("-t", "--threads",          type=int,            metavar="", default=2,        help="Number of sub-processes/threads to use. Default: 2")
        update_group_optional.add_argument("-s", "--specialization",   type=str,            metavar="", default="",       help="Change specialization mode. Can only be used if database was built with some specialization. Options: [sequence,file,assembly,custom]. 'sequence' will use sequence accession as target. 'file' uses the filename as target. 'assembly' will use assembly info from NCBI as target. 'custom' uses the 4th column of the file provided in --seq-info-file as target.")
        update_group_optional.add_argument("--seq-info-mode",          type=str, nargs="*", metavar="", default=["auto"], help="Automatic mode to retrieve tax. info and seq. length. [auto,eutils] or one or more accession2taxid files from NCBI [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot]. auto will either use eutils for less than 50000 input sequences or nucl_gb nucl_wgs. Alternatively a file can be directly provided (see --seq-info-file). Default: auto")
        update_group_optional.add_argument("--seq-info-file",          type=str,            metavar="",                   help="Tab-separated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly id]) [Mutually exclusive --seq-info]")
        update_group_optional.add_argument("--taxdump-file",           type=str, nargs="*", metavar="",                   help="Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)")
        update_group_optional.add_argument("--input-directory",        type=str,            metavar="", default="",       help="Directory containing input files")
        update_group_optional.add_argument("--input-extension",        type=str,            metavar="", default="",       help="Extension of files to use with --input-directory (provide it without * expansion, e.g. '.fna.gz')")
        update_group_optional.add_argument("--update-complete",        action="store_true",                               help="Update adding and removing sequences. Input files should represent the complete updated set of references, not only new sequences.")
        update_group_optional.add_argument("--write-seq-info-file",    action="store_true",                               help="Write sequence information to DB_PREFIX.seqinfo.txt")
        update_group_optional.add_argument("--verbose",                action="store_true",                               help="Verbose output mode")
        update_group_optional.add_argument("--quiet",                  action="store_true",                               help="Quiet output mode")
        update_group_optional.add_argument("--ganon-path",             type=str,            metavar="", default="",       help=argparse.SUPPRESS)
        update_group_optional.add_argument("--n-refs",                 type=int,            metavar="",                   help=argparse.SUPPRESS)
        update_group_optional.add_argument("--n-batches",              type=int,            metavar="",                   help=argparse.SUPPRESS)

        ####################################################################################################

        classify_parser = argparse.ArgumentParser(add_help=False)

        # Required
        classify_group_required = classify_parser.add_argument_group("required arguments")
        classify_group_required.add_argument("-d", "--db-prefix",    type=str, nargs="*", required=True,                                             help="Database input prefix[es]")
        classify_group_required.add_argument("-s", "--single-reads", type=str, nargs="*", required=False, metavar="reads.fq[.gz]",                   help="Multi-fastq[.gz] file[s] to classify")
        classify_group_required.add_argument("-p", "--paired-reads", type=str, nargs="*", required=False, metavar="reads.1.fq[.gz] reads.2.fq[.gz]", help="Multi-fastq[.gz] pairs of file[s] to classify")

        classify_group_cutoff_filter = classify_parser.add_argument_group("cutoff/filter arguments")
        classify_group_cutoff_filter.add_argument("-c", "--rel-cutoff",          type=float, nargs="*", metavar="", default=[0.2],  help="Min. relative percentage of k-mers necessary to consider a match. Generally used to cutoff low similarity matches. Single value or one per database (e.g. 0.5 0.7 1 0.25). 0 for no cutoff. [Mutually exclusive --abs-cutoff] Default: 0.5")
        classify_group_cutoff_filter.add_argument("-e", "--rel-filter",          type=float, nargs="*", metavar="", default=[0.1],  help="Additional relative percentage of k-mers (relative to the best match) to keep a match (applied after cutoff). Single value or one per hierarchy (e.g. 0.1 0 0.25). 1 for no filter. [Mutually exclusive --abs-filter] Default: 0.1")

        classify_group_output = classify_parser.add_argument_group("output arguments")
        classify_group_output.add_argument("-o", "--output-prefix",       type=str,              metavar="", help="Output prefix to print report (.rep). Empty to output to STDOUT")
        classify_group_output.add_argument("--output-lca",                action="store_true",               help="Output an additional file with one lca match for each read (.lca)")
        classify_group_output.add_argument("--output-all",                action="store_true",               help="Output an additional file with all matches. File can be very large (.all)")
        classify_group_output.add_argument("--output-unclassified",       action="store_true",               help="Output an additional file with unclassified read headers (.unc)")
        classify_group_output.add_argument("--output-single",             action="store_true",               help="When using multiple hierarchical levels, output everything in one file instead of one per hierarchy")
        
        classify_group_other = classify_parser.add_argument_group("other arguments")
        classify_group_other.add_argument("-b", "--abs-cutoff",          type=int,   nargs="*", metavar="", help="Max. absolute number of errors allowed to consider a match. Generally used to cutoff low similarity matches. Single value or one per database (e.g. 3 3 4 0). -1 for no cutoff. [Mutually exclusive --rel-cutoff]")
        classify_group_other.add_argument("-a", "--abs-filter",          type=int,   nargs="*", metavar="", help="Additional absolute number of errors (relative to the best match) to keep a match (applied after cutoff). Single value or one per hierarchy (e.g. 0 2 1). -1 for no filter. [Mutually exclusive --rel-filter]")
        classify_group_other.add_argument("-l", "--hierarchy-labels",    type=str,   nargs="*", metavar="", help="Hierarchy definition, one for each database input. Can also be a string, but input will be sorted to define order (e.g. 1 1 2 3). Default: H1")
        classify_group_other.add_argument("-f", "--offset",              type=int,              metavar="", help="Number of k-mers to skip during classification. Can speed up analysis but may reduce recall. (e.g. 1 = all k-mers, 3 = every 3rd k-mer). Default: 1")
        classify_group_other.add_argument("-t", "--threads",             type=int,              metavar="", help="Number of sub-processes/threads to use. Default: 3")
        classify_group_other.add_argument("-r", "--ranks",               type=str,   nargs="*", metavar="", help="Ranks to show in the report (.tre). 'all' for all identified ranks. empty for default ranks: superkingdom phylum class order family genus species assembly. This file can be re-generated with the ganon report command.")
        classify_group_other.add_argument("--verbose",                   action="store_true",               help="Verbose output mode")
        classify_group_other.add_argument("--quiet",                     action="store_true",               help="Quiet output mode")
        classify_group_other.add_argument("--ganon-path",                type=str, default="",  metavar="", help=argparse.SUPPRESS) 
        classify_group_other.add_argument("--n-reads",                   type=int,              metavar="", help=argparse.SUPPRESS)
        classify_group_other.add_argument("--n-batches",                 type=int,              metavar="", help=argparse.SUPPRESS)
        classify_group_other.add_argument("--hibf" ,                     action="store_true",               help=argparse.SUPPRESS)

        ####################################################################################################

        report_parser = argparse.ArgumentParser(add_help=False)

        # Required
        report_group_required = report_parser.add_argument_group("required arguments")
        report_group_required.add_argument("-i", "--rep-files",     type=str, nargs="*", required=False, help="One or more *.rep files from ganon classify")
        report_group_required.add_argument("-o", "--output-prefix", type=str,            required=True,  help="Output prefix for report file '{output_prefix}.tre'. In case of multiple files, the base input filename will be appended at the end of the output file '{output_prefix + FILENAME}.tre'")

        # Defaults
        report_group_optional = report_parser.add_argument_group("optional arguments")
        report_group_optional.add_argument("-d", "--db-prefix",      type=str, nargs="*", metavar="", default=[],      help="Database prefix[es] used for classification (in any order). Only '.tax' file is required. If not provided, new taxonomy will be downloaded")
        report_group_optional.add_argument("-f", "--output-format",  type=str,            metavar="", default="tsv",   help="Output format [text, tsv, csv]. text outputs a tabulated formatted text file for better visualization. Default: tsv")
        report_group_optional.add_argument("-e", "--report-type",    type=str,            metavar="", default="reads", help="Type of report to generate [reads, matches]. Default: reads")
        report_group_optional.add_argument("-r", "--ranks",          type=str, nargs="*", metavar="", default=[],      help="Ranks to report ['', 'all', custom list] 'all' for all possible ranks. empty for default ranks (superkingdom phylum class order family genus species assembly). Default: """)
        report_group_optional.add_argument("-s", "--sort",           type=str,            metavar="", default="",      help="Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage (with --ranks all)")
        report_group_optional.add_argument("-y", "--split-hierarchy",action="store_true",                              help="Split output reports by hierarchy (from ganon classify --hierarchy-labels). If activated, the output files will be named as '{output_prefix}.{hierarchy}.tre'")
        report_group_optional.add_argument("-p", "--skip-hierarchy", type=str, nargs="*", metavar="", default=[],      help="One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels)")
        report_group_optional.add_argument("-k", "--keep-hierarchy", type=str, nargs="*", metavar="", default=[],      help="One or more hierarchies to keep in the report (from ganon classify --hierarchy-labels)")
        report_group_optional.add_argument("--taxdump-file",         type=str, nargs="*", metavar="", default=[],      help="Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp]) file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)")
        report_group_optional.add_argument("--input-directory",      type=str,            metavar="", default="",      help="Directory containing input files")
        report_group_optional.add_argument("--input-extension",      type=str,            metavar="", default="",      help="Extension of files to use with --input-directory (provide it without * expansion, e.g. '.rep')")
        report_group_optional.add_argument("--verbose",              action="store_true",             default=False,   help="Verbose output mode")
        report_group_optional.add_argument("--quiet",                action="store_true",             default=False,   help="Quiet output mode")

        ####################################################################################################

        table_parser = argparse.ArgumentParser(add_help=False)

        # Required
        table_group_required = table_parser.add_argument_group("required arguments")
        table_group_required.add_argument("-i", "--tre-files",   type=str, nargs="*", required=False, help="Report files (.tre) from ganon classify/report to make the table")
        table_group_required.add_argument("-o", "--output-file", type=str,            required=True,  help="Output filename for the table")

        # Defaults
        table_group_optional = table_parser.add_argument_group("optional arguments")
        table_group_optional.add_argument("-l", "--output-value",             type=str,   metavar="", default="counts",     help="Output value on the table [percentage, counts]. percentage values are reported between [0-1]. Default: counts")
        table_group_optional.add_argument("-f", "--output-format",            type=str,   metavar="", default="tsv",        help="Output format [tsv, csv]. Default: tsv")
        table_group_optional.add_argument("-t", "--top-sample",               type=int,   metavar="", default=0,            help="Top hits of each sample individually")
        table_group_optional.add_argument("-a", "--top-all",                  type=int,   metavar="", default=0,            help="Top hits of all samples (ranked by percentage)")
        table_group_optional.add_argument("-m", "--min-frequency",            type=float, metavar="", default=0,            help="Minimum number/percentage of files containing an taxa to keep the taxa [values between 0-1 for percentage, >1 specific number]")
        table_group_optional.add_argument("-r", "--rank",                     type=str,   metavar="", default=None,         help="Define specific rank to report. Empty will report all ranks.")
        table_group_optional.add_argument("-n", "--no-root",                  action="store_true",    default=False,        help="Do not report root node entry and lineage. Direct and shared matches to root will be accounted as unclassified")
        table_group_optional.add_argument("--header",                         type=str,   metavar="", default="name",       help="Header information [name, taxid, lineage]. Default: name")
        table_group_optional.add_argument("--unclassified-label",             type=str,   metavar="", default=None,         help="Add column with unclassified count/percentage with the chosen label. May be the same as --filtered-label (e.g. unassigned)")
        table_group_optional.add_argument("--filtered-label",                 type=str,   metavar="", default=None,         help="Add column with filtered count/percentage with the chosen label. May be the same as --unclassified-label (e.g. unassigned)")
        table_group_optional.add_argument("--skip-zeros",                     action="store_true",    default=False,        help="Do not print lines with only zero count/percentage")
        table_group_optional.add_argument("--transpose",                      action="store_true",    default=False,        help="Transpose output table (taxa as cols and files as rows)")
        table_group_optional.add_argument("--input-directory",                type=str,  metavar="",  default="",           help="Directory containing input files")
        table_group_optional.add_argument("--input-extension",                type=str,  metavar="",  default="",           help="Extension of files to use with --input-directory (provide it without * expansion, e.g. '.tre')")
        table_group_optional.add_argument("--verbose",                        action="store_true",    default=False,        help="Verbose output mode")
        table_group_optional.add_argument("--quiet",                          action="store_true",    default=False,        help="Quiet output mode")

        ####################################################################################################

        filter_parser = argparse.ArgumentParser(add_help=False)
        filter_arguments = filter_parser.add_argument_group("filter arguments")
        filter_arguments.add_argument("--min-count",      type=float,          metavar="", default=0,  help="Minimum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1 specific number]")
        filter_arguments.add_argument("--max-count",      type=float,          metavar="", default=0,  help="Maximum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1 specific number]")
        filter_arguments.add_argument("--names",          type=str, nargs="*", metavar="", default=[], help="Show only entries matching exact names of the provided list")
        filter_arguments.add_argument("--names-with",     type=str, nargs="*", metavar="", default=[], help="Show entries containing full or partial names of the provided list")
        filter_arguments.add_argument("--taxids",         type=str, nargs="*", metavar="", default=[], help="One or more taxids to report (including children taxa)")

        formatter_class = lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=120)
        subparsers = parser.add_subparsers()

        build = subparsers.add_parser("build", help="Download and build ganon default databases (refseq/genbank)", parents=[build_parser, build_default_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        #build._action_groups.reverse()  # required first
        build.set_defaults(which="build")

        update = subparsers.add_parser("update", help="Update ganon default databases", parents=[update_parser], formatter_class=formatter_class)
        update._action_groups.reverse()  # required first
        update.set_defaults(which="update")

        classify = subparsers.add_parser("classify", help="Classify reads", parents=[classify_parser], formatter_class=formatter_class)
        classify.set_defaults(which="classify")

        report = subparsers.add_parser("report", help="Generate reports", parents=[filter_parser,report_parser], formatter_class=formatter_class)
        report._action_groups.reverse()  # required first
        report.set_defaults(which="report")

        table = subparsers.add_parser("table", help="Generate table from reports", parents=[filter_parser,table_parser], formatter_class=formatter_class)
        table._action_groups.reverse()  # required first
        table.set_defaults(which="table")

        build_custom = subparsers.add_parser("build-custom", help="Build custom ganon databases", parents=[build_default_parser, build_custom_parser], formatter_class=formatter_class)
        build_custom.set_defaults(which="build-custom")

        update_custom = subparsers.add_parser("update-custom", help="Update custom ganon databases", parents=[update_custom_parser], formatter_class=formatter_class)
        update_custom.set_defaults(which="update-custom")

        # Passing arguments internally from call main(which, **kwargs)
        if which is not None:
            # Set which as the first parameter (mandatory)
            list_kwargs = [which]
            for arg, value in kwargs.items():
                # convert all others to argparse format
                # (eg: input_files to --input-files)
                arg_formatted = "--" + arg.replace("_", "-")
                if isinstance(value, list):  # unpack if list
                    list_kwargs.append(arg_formatted)
                    list_kwargs.extend(value)
                elif type(value) == bool and value is True:  # add only arg if boolean flag activated
                    list_kwargs.append(arg_formatted)
                elif value:
                    list_kwargs.append(arg_formatted)
                    list_kwargs.append(str(value))
            # Parse from list saving arguments to this class
            parser.parse_args(list_kwargs, namespace=self)
        else:
            # parse from default CLI sys.argv saving arguments to this class
            parser.parse_args(namespace=self)
            if len(sys.argv) == 1:
                parser.print_help()
                self.empty = True

    def __repr__(self):
        args = ["{}={}".format(k, repr(v)) for (k,v) in vars(self).items()]
        return "Config({})".format(", ".join(args))

    def validate(self):

        if self.empty is True:
            print_log("Please provide one or more arguments")
            return False

        if self.verbose is True:
            self.quiet = False
        elif self.quiet is True:
            self.verbose = False

        if self.which == "build":
            pass

        elif self.which == "build-custom":

            if self.input_file and self.input:
                print_log("--input-file is mutually exclusive with --input")
                return False

            if self.level == "custom" and not self.input_file:
                print_log("--level custom requires --input-file")
                return False

            if self.level and self.level not in self.choices_level and self.taxonomy == "none":
                print_log("--taxonomy is required for --level " + self.level)
                return False

            if self.taxonomy=="ncbi":
                for entry in self.ncbi_sequence_info:
                    if entry not in self.choices_ncbi_sequence_info and not check_file(entry):
                        print_log("Invalid --get-sequence-info. Should be a valid file or: " + " ".join(self.choices_ncbi_sequence_info))
                        return False

                for entry in self.ncbi_file_info:
                    if entry not in self.choices_ncbi_file_info and not check_file(entry):
                        print_log("Invalid --get-file-info. Should be a valid file or: " + " ".join(self.choices_ncbi_file_info))
                        return False


        elif self.which == "update" or self.which == "update-custom":
            if not check_file(self.db_prefix+".ibf"):
                return False

            if not all([sim in retrieve_info_options for sim in self.retrieve_info]):
                print_log("Invalid --seq-info-mode. Options: " + " ".join(seq_info_mode_options))
                return False

        elif self.which == "classify":
            if not all([check_file(prefix + ".ibf") for prefix in self.db_prefix]):
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

            if len_single_reads + len_paired_reads==0:
                print_log("No valid input files to classify")
                return False

        elif self.which=="report":

            if self.skip_hierarchy and self.keep_hierarchy:
                print_log("--skip-hierarchy and --keep-hierarchy are mutually exclusive")
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

        elif self.which == "table":
            if self.min_frequency < 0:
                print_log("Invalid value for --min-frequency (>0)")
                return False

        return True

    def set_paths(self):
        missing_path = False
        if self.which in ["build", "build-custom", "update", "update-custom"]:
            self.ganon_path = self.ganon_path + "/" if self.ganon_path else ""

            # if path is given, look for binaries only there
            ganon_build_paths = [self.ganon_path, self.ganon_path+"build/"] if self.ganon_path else [None, "build/"]
            for p in ganon_build_paths:
                self.path_exec["build"] = shutil.which("ganon-build", path=p)
                if self.path_exec["build"] is not None: break
            if self.path_exec["build"] is None:
                print_log("ganon-build binary was not found. Please inform a specific path with --ganon-path")
                missing_path = True

            ganon_get_seq_info_paths = [self.ganon_path, self.ganon_path+"scripts/", self.ganon_path+"../scripts/"] if self.ganon_path else [None, "scripts/"]
            for p in ganon_get_seq_info_paths:
                self.path_exec["get_seq_info"] = shutil.which("ganon-get-seq-info.sh", path=p)
                if self.path_exec["get_seq_info"] is not None: break
            if self.path_exec["get_seq_info"] is None:
                print_log("ganon-get-seq-info.sh script was not found. Please inform a specific path with --ganon-path")
                missing_path = True

            ganon_genome_updater_paths = [self.ganon_path, self.ganon_path+"libs/genome_updater/", self.ganon_path+"../libs/genome_updater/"] if self.ganon_path else [None, "libs/genome_updater/"]
            for p in ganon_genome_updater_paths:
                self.path_exec["genome_updater"] = shutil.which("genome_updater.sh", path=p)
                if self.path_exec["genome_updater"] is not None: break
            if self.path_exec["genome_updater"] is None:
                print_log("genome_updater.sh was not found. Please inform a specific path with --ganon-path")
                missing_path = True

        elif self.which in ["classify"]:
            self.ganon_path = self.ganon_path + "/" if self.ganon_path else ""

            ganon_classify_paths = [self.ganon_path, self.ganon_path+"build/"] if self.ganon_path else [None, "build/"]
            for p in ganon_classify_paths:
                self.path_exec["classify"] = shutil.which("ganon-classify", path=p)
                if self.path_exec["classify"] is not None: break
            if self.path_exec["classify"] is None:
                print_log("ganon-classify binary was not found. Please inform a specific path with --ganon-path")
                missing_path = True

        return True if not missing_path else False


def int_or_float(minval=None, maxval=None):
    def checker(val):
        try:
            val = int(val)
        except ValueError:
            val = float(val)
        if minval is not None and val < minval:
            raise argparse.ArgumentTypeError('%r must be >= %r' % (val, minval))
        if maxval is not None and val > maxval:
            raise argparse.ArgumentTypeError('%r must be <= %r' % (val, maxval))
        return val
    return checker

def file_exists(file):
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError("{0} does not exist".format(file))
    return file