#!/usr/bin/env python3

import argparse
import sys
import shutil
from ganon.util import *


class Config:

    version = "2.1.1"
    path_exec = {"build": "", "classify": "", "get_seq_info": "", "genome_updater": ""}
    empty = False

    choices_taxonomy = ["ncbi", "gtdb", "skip"]  # get from multitax
    choices_og = ["archaea", "bacteria", "fungi", "human", "invertebrate", "metagenomes",
                  "other", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"]
    choices_db_source = ["refseq", "genbank"]
    choices_level = ["assembly", "custom"]
    choices_input_target = ["file", "sequence"]
    choices_ncbi_sequence_info = ["eutils", "nucl_gb", "nucl_wgs", "nucl_est", "nucl_gss", "pdb",
                                  "prot", "dead_nucl", "dead_wgs", "dead_prot"]
    choices_ncbi_file_info = ["refseq", "genbank", "refseq_historical", "genbank_historical"]
    choices_default_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "assembly"]
    choices_report_type = ["abundance", "reads", "matches", "dist", "corr"]
    choices_multiple_matches = ["em", "lca", "skip"]
    choices_report_output = ["text", "tsv", "csv", "bioboxes"]
    choices_mode = ["avg", "smaller", "smallest", "faster", "fastest"]
    choices_filter_type = ["hibf", "ibf"]

    def __init__(self, which: str=None, **kwargs):

        parser = argparse.ArgumentParser(prog="ganon",
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description=logo(self.version),
                                         conflict_handler="resolve")
        parser.add_argument("-v", "--version", action="version", version="version: %(prog)s " + self.version, help="Show program's version number and exit.")

        ####################################################################################################

        build_default_parser = argparse.ArgumentParser(add_help=False)

        build_default_required_args = build_default_parser.add_argument_group("required arguments")
        build_default_required_args.add_argument("-d", "--db-prefix", type=str, required=True, help="Database output prefix")

        build_default_important_args = build_default_parser.add_argument_group("important arguments")
        build_default_important_args.add_argument("-x", "--taxonomy", type=str,                    metavar="", default="ncbi", help="Set taxonomy to enable taxonomic classification, lca and reports [" + ", ".join(self.choices_taxonomy) + "]", choices=self.choices_taxonomy)
        build_default_important_args.add_argument("-t", "--threads",  type=unsigned_int(minval=1), metavar="", default=1,      help="")

        build_default_advanced_args = build_default_parser.add_argument_group("advanced arguments")
        build_default_advanced_args.add_argument("-p", "--max-fp",         type=int_or_float(minval=0, maxval=1), metavar="", default=None,   help="Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with --filter-type hibf or 0.05 with --filter-type ibf.")
        build_default_advanced_args.add_argument("-k", "--kmer-size",      type=unsigned_int(minval=1),           metavar="", default=19,     help="The k-mer size to split sequences.")
        build_default_advanced_args.add_argument("-w", "--window-size",    type=unsigned_int(minval=1),           metavar="", default=31,     help="The window-size to build filter with minimizers.")
        build_default_advanced_args.add_argument("-s", "--hash-functions", type=unsigned_int(minval=0, maxval=5), metavar="", default=4,      help="The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0 will try to set optimal value.", choices=range(6))
        build_default_advanced_args.add_argument("-f", "--filter-size",    type=unsigned_float(),                 metavar="", default=0,      help="Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-type ibf.")
        build_default_advanced_args.add_argument("-j", "--mode",           type=str,                              metavar="", default="avg",  help="Create smaller or faster filters at the cost of classification speed or database size, respectively [" + ", ".join(self.choices_mode) + "]. If --filter-size is used, smaller/smallest refers to the false positive rate. By default, an average value is calculated to balance classification speed and database size. Only valid for --filter-type ibf.", choices=self.choices_mode)
        build_default_advanced_args.add_argument("-y", "--min-length",     type=unsigned_int(minval=0),           metavar="", default=0,      help="Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-type ibf.")
        build_default_advanced_args.add_argument("-v", "--filter-type",    type=str,                              metavar="", default="hibf", help="Variant of bloom filter to use [" + ", ".join(self.choices_filter_type) + "]. hibf requires raptor >= v3.0.1 installed or binary path set with --raptor-path. --mode, --filter-size and --min-length will be ignored with hibf. hibf will set --max-fp 0.001 as default.", choices=self.choices_filter_type)

        ####################################################################################################

        build_parser = argparse.ArgumentParser(add_help=False)

        build_required_args = build_parser.add_argument_group("required arguments")
        build_required_args.add_argument("-g", "--organism-group", type=str, nargs="*", metavar="", help="One or more organism groups to download [" + ", ".join(self.choices_og) + "]. Mutually exclusive --taxid", choices=self.choices_og)
        build_required_args.add_argument("-a", "--taxid",          type=str, nargs="*", metavar="", help="One or more taxonomic identifiers to download. e.g. 562 (-x ncbi) or 's__Escherichia coli' (-x gtdb). Mutually exclusive --organism-group")  

        build_database_args = build_parser.add_argument_group("database arguments")
        build_database_args.add_argument("-l", "--level",          type=str, default="species", metavar="", help="Highest level to build the database. Options: any available taxonomic rank [species, genus, ...], 'leaves' for taxonomic leaves or 'assembly' for a assembly/strain based analysis")
        
        build_download_args = build_parser.add_argument_group("download arguments")
        build_download_args.add_argument("-b", "--source",            type=str, nargs="*",         default=["refseq"], metavar="", help="Source to download [" + ", ".join(self.choices_db_source) + "]", choices=self.choices_db_source)
        build_download_args.add_argument("-o", "--top",               type=unsigned_int(minval=0), default=0,          metavar="", help="Download limited assemblies for each taxa. 0 for all.")
        build_download_args.add_argument("-c", "--complete-genomes",        action="store_true",                                         help="Download only sub-set of complete genomes")
        build_download_args.add_argument("-r", "--representative-genomes",  action="store_true",                                         help="Download only sub-set of representative genomes")
        build_download_args.add_argument("-u", "--genome-updater",    type=str,                                        metavar="", help="Additional genome_updater parameters (https://github.com/pirovc/genome_updater)")
        build_download_args.add_argument("-m", "--taxonomy-files",    type=file_exists, nargs="*", metavar="",                     help="Specific files for taxonomy - otherwise files will be downloaded")
        build_download_args.add_argument("-z", "--genome-size-files", type=file_exists, nargs="*", metavar="",                     help="Specific files for genome size estimation - otherwise files will be downloaded")
        build_download_args.add_argument("--skip-genome-size",        action="store_true",  help="Do not attempt to get genome sizes. Activate this option when using sequences not representing full genomes.")

        ####################################################################################################

        build_custom_parser = argparse.ArgumentParser(add_help=False)

        build_custom_required_args = build_custom_parser.add_argument_group("required arguments")
        build_custom_required_args.add_argument("-i", "--input",           type=str,          nargs="*",        metavar="", help="Input file(s) and/or folder(s). Mutually exclusive --input-file.")
        build_custom_required_args.add_argument("-e", "--input-extension", type=str,          default="fna.gz", metavar="", help="Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).")
        build_custom_required_args.add_argument("-c", "--input-recursive", action="store_true",                             help="Look for files recursively in folder(s) provided with --input")

        build_custom_args = build_custom_parser.add_argument_group("custom arguments")
        build_custom_args.add_argument("-n", "--input-file",        type=file_exists,                 metavar="", help="Tab-separated file with all necessary file/sequence information. Fields: file [<tab> target <tab> node <tab> specialization <tab> specialization name]. For details: https://pirovc.github.io/ganon/custom_databases/. Mutually exclusive --input")
        build_custom_args.add_argument("-a", "--input-target",      type=str,         default="file", metavar="", help="Target to use [file, sequence]. Parse input by file or by sequence. Using 'file' is recommended and will speed-up the building process", choices=self.choices_input_target)
        build_custom_args.add_argument("-l", "--level",             type=str,                         metavar="", help="Max. level to build the database. By default, --level is the --input-target. Options: any available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy). Further specialization options [" + ", ".join(self.choices_level) + "]. assembly will retrieve and use the assembly accession and name. custom requires and uses the specialization field in the --input-file.")
        build_custom_args.add_argument("-m", "--taxonomy-files",    type=file_exists, nargs="*",      metavar="", help="Specific files for taxonomy - otherwise files will be downloaded")
        build_custom_args.add_argument("-z", "--genome-size-files", type=file_exists, nargs="*",      metavar="", help="Specific files for genome size estimation - otherwise files will be downloaded")
        build_custom_args.add_argument("--skip-genome-size",        action="store_true",  help="Do not attempt to get genome sizes. Activate this option when using sequences not representing full genomes.")

        ncbi_args = build_custom_parser.add_argument_group("ncbi arguments")
        ncbi_args.add_argument("-r", "--ncbi-sequence-info", type=str, nargs="*", default=[],                               metavar="", help="Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information. [" + ", ".join(self.choices_ncbi_sequence_info) + " or one or more accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/]. By default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise.")
        ncbi_args.add_argument("-q", "--ncbi-file-info",     type=str, nargs="*", default=self.choices_ncbi_file_info[0:2], metavar="", help="Downloads assembly_summary files to extract target information. [" + ", ".join(self.choices_ncbi_file_info) + " or one or more assembly_summary files from https://ftp.ncbi.nlm.nih.gov/genomes/]")

        ####################################################################################################

        update_parser = argparse.ArgumentParser(add_help=False)

        # Required
        update_group_required = update_parser.add_argument_group("required arguments")
        update_group_required.add_argument("-d", "--db-prefix", type=str, required=True, help="Existing database input prefix")

        update_default_important_args = update_parser.add_argument_group("important arguments")
        update_default_important_args.add_argument("-o", "--output-db-prefix", type=str,            metavar="",            help="Output database prefix. By default will be the same as --db-prefix and overwrite files")
        update_default_important_args.add_argument("-t", "--threads",  type=unsigned_int(minval=1), metavar="", default=1, help="")

        ####################################################################################################

        build_update_parser = argparse.ArgumentParser(add_help=False)

        build_update_other_args = build_update_parser.add_argument_group("optional arguments")
        build_update_other_args.add_argument("--restart",         action="store_true", help="Restart build/update from scratch, do not try to resume from the latest possible step. {db_prefix}_files/ will be deleted if present.")
        build_update_other_args.add_argument("--verbose",         action="store_true", help="Verbose output mode")
        build_update_other_args.add_argument("--quiet",           action="store_true", help="Quiet output mode")
        build_update_other_args.add_argument("--keep-files",      action="store_true", help=argparse.SUPPRESS)
        build_update_other_args.add_argument("--write-info-file", action="store_true", help="Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for further attempts.")
        build_update_other_args.add_argument("--ganon-path",  type=str,                    metavar="", default="", help=argparse.SUPPRESS)
        build_update_other_args.add_argument("--raptor-path", type=str,                    metavar="", default="", help=argparse.SUPPRESS)
        build_update_other_args.add_argument("--n-refs",      type=unsigned_int(minval=1), metavar="",             help=argparse.SUPPRESS)
        build_update_other_args.add_argument("--n-batches",   type=unsigned_int(minval=1), metavar="",             help=argparse.SUPPRESS)
        build_update_other_args.add_argument("--ncbi-url",    type=str,                    metavar="", default="https://ftp.ncbi.nlm.nih.gov/", help=argparse.SUPPRESS)
        build_update_other_args.add_argument("--gtdb-url",    type=str,                    metavar="", default="https://data.gtdb.ecogenomic.org/releases/latest/", help=argparse.SUPPRESS)

        ####################################################################################################

        classify_parser = argparse.ArgumentParser(add_help=False)

        # Required
        classify_group_required = classify_parser.add_argument_group("required arguments")
        classify_group_required.add_argument("-d", "--db-prefix",    type=str, nargs="*", required=True,                                             help="Database input prefix[es]")
        classify_group_required.add_argument("-s", "--single-reads", type=str, nargs="*", required=False, metavar="reads.fq[.gz]",                   help="Multi-fastq[.gz] file[s] to classify")
        classify_group_required.add_argument("-p", "--paired-reads", type=str, nargs="*", required=False, metavar="reads.1.fq[.gz] reads.2.fq[.gz]", help="Multi-fastq[.gz] pairs of file[s] to classify")

        classify_group_cutoff_filter = classify_parser.add_argument_group("cutoff/filter arguments")
        classify_group_cutoff_filter.add_argument("-c", "--rel-cutoff", type=int_or_float(minval=0, maxval=1), nargs="*", metavar="", default=[0.75], help="Min. percentage of a read (set of k-mers) shared with a reference necessary to consider a match. Generally used to remove low similarity matches. Single value or one per database (e.g. 0.7 1 0.25). 0 for no cutoff")
        classify_group_cutoff_filter.add_argument("-e", "--rel-filter", type=int_or_float(minval=0, maxval=1), nargs="*", metavar="", default=[0.1],  help="Additional relative percentage of matches (relative to the best match) to keep. Generally used to keep top matches above cutoff. Single value or one per hierarchy (e.g. 0.1 0). 1 for no filter")
      
        classify_group_postrep = classify_parser.add_argument_group("post-processing/report arguments")
        classify_group_postrep.add_argument("-m", "--multiple-matches", type=str,                    metavar="", default="em",        help="Method to solve reads with multiple matches  [" + ", ".join(self.choices_multiple_matches) + "]. em -> expectation maximization algorithm based on unique matches. lca -> lowest common ancestor based on taxonomy. The EM algorithm can be executed later with 'ganon reassign' using the .all file (--output-all).", choices=self.choices_multiple_matches)
        classify_group_postrep.add_argument("--ranks",                  type=str,  nargs="*",        metavar="", default=[],          help="Ranks to report taxonomic abundances (.tre). empty will report default ranks [" + ", ".join(self.choices_default_ranks) + "].")
        classify_group_postrep.add_argument("--min-count",              type=int_or_float(minval=0), metavar="", default=0.00005,     help="Minimum percentage/counts to report an taxa (.tre) [use values between 0-1 for percentage, >1 for counts]")
        classify_group_postrep.add_argument("--report-type",            type=str,                    metavar="", default="abundance", help="Type of report (.tre) [" + ", ".join(self.choices_report_type) + "]. More info in 'ganon report'.", choices=self.choices_report_type)
        classify_group_postrep.add_argument("--skip-report",            action="store_true",                                          help="Disable tree-like report (.tre) at the end of classification. Can be done later with 'ganon report'.")

        classify_group_output = classify_parser.add_argument_group("output arguments")
        classify_group_output.add_argument("-o", "--output-prefix",       type=str,              metavar="", help="Output prefix for output (.rep) and tree-like report (.tre). Empty to output to STDOUT (only .rep)")
        classify_group_output.add_argument("--output-one",                action="store_true",               help="Output a file with one match for each read (.one) either an unique match or a result from the EM or a LCA algorithm (--multiple-matches)")
        classify_group_output.add_argument("--output-all",                action="store_true",               help="Output a file with all unique and multiple matches (.all)")
        classify_group_output.add_argument("--output-unclassified",       action="store_true",               help="Output a file with unclassified read headers (.unc)")
        classify_group_output.add_argument("--output-single",             action="store_true",               help="When using multiple hierarchical levels, output everything in one file instead of one per hierarchy")

        classify_group_other = classify_parser.add_argument_group("other arguments")
        classify_group_other.add_argument("-t", "--threads",             type=unsigned_int(minval=1), metavar="", default=1,  help="Number of sub-processes/threads to use")
        classify_group_other.add_argument("-b", "--binning",             action="store_true",                                 help="Optimized parameters for binning (--rel-cutoff 0.25 --rel-filter 0 --min-count 0 --report-type reads). Will report sequence abundances (.tre) instead of tax. abundance.")
        classify_group_other.add_argument("-f", "--fpr-query",  type=int_or_float(minval=0, maxval=1), nargs="*", metavar="", default=[1e-5], help="Max. false positive of a query to accept a match. Applied after --rel-cutoff and --rel-filter. Generally used to remove false positives matches querying a database build with large --max-fp. Single value or one per hierarchy (e.g. 0.1 0). 1 for no filter")
        classify_group_other.add_argument("-l", "--hierarchy-labels",    type=str,         nargs="*", metavar="",             help="Hierarchy definition of --db-prefix files to be classified. Can also be a string, but input will be sorted to define order (e.g. 1 1 2 3). The default value reported without hierarchy is 'H1'")
        classify_group_other.add_argument("--verbose",                   action="store_true",                                 help="Verbose output mode")
        classify_group_other.add_argument("--quiet",                     action="store_true",                                 help="Quiet output mode")
        classify_group_other.add_argument("--hibf",                      action="store_true",                     help=argparse.SUPPRESS)
        classify_group_other.add_argument("--ganon-path",                type=str, default="",  metavar="",       help=argparse.SUPPRESS) 
        classify_group_other.add_argument("--n-reads",                   type=unsigned_int(minval=1), metavar="", help=argparse.SUPPRESS)
        classify_group_other.add_argument("--n-batches",                 type=unsigned_int(minval=1), metavar="", help=argparse.SUPPRESS)

        ####################################################################################################

        reassign_parser = argparse.ArgumentParser(add_help=False)

        # Required
        reassign_group_required = reassign_parser.add_argument_group("required arguments")
        reassign_group_required.add_argument("-i", "--input-prefix",  type=str, required=True, metavar="", help="Input prefix to find files from ganon classify (.all and optionally .rep)")
        reassign_group_required.add_argument("-o", "--output-prefix", type=str, required=True,             help="Output prefix for reassigned file (.one and optionally .rep). In case of multiple files, the base input filename will be appended at the end of the output file 'output_prefix + FILENAME.out'")
   
        reassign_em = reassign_parser.add_argument_group("EM arguments")
        reassign_em.add_argument("-e", "--max-iter",  type=unsigned_int(minval=0), metavar="", default=10, help="Max. number of iterations for the EM algorithm. If 0, will run until convergence (check --threshold)")
        reassign_em.add_argument("-s", "--threshold", type=int_or_float(minval=0), metavar="", default=0,  help="Convergence threshold limit to stop the EM algorithm.")

        reassign_group_other = reassign_parser.add_argument_group("other arguments")
        reassign_group_other.add_argument("--remove-all", action="store_true", help="Remove input file (.all) after processing.")
        reassign_group_other.add_argument("--skip-one",   action="store_true", help="Do not write output file (.one) after processing.")
        reassign_group_other.add_argument("--verbose",    action="store_true", help="Verbose output mode")
        reassign_group_other.add_argument("--quiet",      action="store_true", help="Quiet output mode")

        ####################################################################################################

        report_parser = argparse.ArgumentParser(add_help=False)

        report_group_required = report_parser.add_argument_group("required arguments")
        report_group_required.add_argument("-i", "--input",           type=str, required=True, nargs="*", metavar="", help="Input file(s) and/or folder(s). '.rep' file(s) from ganon classify.")
        report_group_required.add_argument("-e", "--input-extension", type=str, default="rep",            help="Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).")
        report_group_required.add_argument("-o", "--output-prefix",   type=str, required=True,            help="Output prefix for report file 'output_prefix.tre'. In case of multiple files, the base input filename will be appended at the end of the output file 'output_prefix + FILENAME.tre'")

        report_group_dbtax = report_parser.add_argument_group("db/tax arguments")
        report_group_dbtax.add_argument("-d", "--db-prefix",         type=str,         nargs="*", metavar="", default=[],     help="Database prefix(es) used for classification. Only '.tax' file(s) are required. If not provided, new taxonomy will be downloaded. Mutually exclusive with --taxonomy.")
        report_group_dbtax.add_argument("-x", "--taxonomy",          type=str,                    metavar="", default="ncbi", help="Taxonomy database to use [" + ", ".join(self.choices_taxonomy) + "]. Mutually exclusive with --db-prefix.", choices=self.choices_taxonomy)
        report_group_dbtax.add_argument("-m", "--taxonomy-files",    type=file_exists, nargs="*", metavar="",                 help="Specific files for taxonomy - otherwise files will be downloaded")
        report_group_dbtax.add_argument("-z", "--genome-size-files", type=file_exists, nargs="*", metavar="",                 help="Specific files for genome size estimation - otherwise files will be downloaded")
        report_group_dbtax.add_argument("--skip-genome-size",        action="store_true",  help="Do not attempt to get genome sizes. Valid only without --db-prefix. Activate this option when using sequences not representing full genomes.")

        report_group_output = report_parser.add_argument_group("output arguments")
        report_group_output.add_argument("-f", "--output-format",  type=str,            metavar="", default="tsv",       help="Output format [" + ", ".join(self.choices_report_output) + "]. text outputs a tabulated formatted text file for better visualization. bioboxes is the the CAMI challenge profiling format (only percentage/abundances are reported).", choices=self.choices_report_output)
        report_group_output.add_argument("-t", "--report-type",    type=str,            metavar="", default="abundance", help="Type of report [" + ", ".join(self.choices_report_type) + "]. 'abundance' -> tax. abundance (re-distribute read counts and correct by genome size), 'reads' -> sequence abundance, 'matches' -> report all unique and shared matches, 'dist' -> like reads with re-distribution of shared read counts only, 'corr' -> like abundance without re-distribution of shared read counts", choices=self.choices_report_type)
        report_group_output.add_argument("-r", "--ranks",          type=str, nargs="*", metavar="", default=[],          help="Ranks to report ['', 'all', custom list]. 'all' for all possible ranks. empty for default ranks [" + ", ".join(self.choices_default_ranks) + "].")
        report_group_output.add_argument("-s", "--sort",           type=str,            metavar="", default="",          help="Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage (with --ranks all)")
        report_group_output.add_argument("-a", "--no-orphan",       action="store_true",                                 help="Omit orphan nodes from the final report. Otherwise, orphan nodes (= nodes not found in the db/tax) are reported as 'na' with root as direct parent.")
        report_group_output.add_argument("-y", "--split-hierarchy", action="store_true",                                 help="Split output reports by hierarchy (from ganon classify --hierarchy-labels). If activated, the output files will be named as '{output_prefix}.{hierarchy}.tre'")
        report_group_output.add_argument("-p", "--skip-hierarchy", type=str,                              nargs="*", metavar="", default=[],          help="One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels)")
        report_group_output.add_argument("-k", "--keep-hierarchy", type=str,                              nargs="*", metavar="", default=[],          help="One or more hierarchies to keep in the report (from ganon classify --hierarchy-labels)")
        report_group_output.add_argument("-c", "--top-percentile", type=int_or_float(minval=0, maxval=0.999999),     metavar="", default=0,           help="Top percentile filter, based on percentage/relative abundance. Applied only at default ranks [" + ", ".join(self.choices_default_ranks) + "]")

        report_group_optional = report_parser.add_argument_group("optional arguments")
        report_group_optional.add_argument("--verbose", action="store_true", default=False, help="Verbose output mode")
        report_group_optional.add_argument("--quiet",   action="store_true", default=False, help="Quiet output mode")
        report_group_optional.add_argument("--ncbi-url",             type=str,                              metavar="", default="https://ftp.ncbi.nlm.nih.gov/", help=argparse.SUPPRESS)
        report_group_optional.add_argument("--gtdb-url",             type=str,                              metavar="", default="https://data.gtdb.ecogenomic.org/releases/latest/", help=argparse.SUPPRESS)

        ####################################################################################################

        table_parser = argparse.ArgumentParser(add_help=False)

        table_group_required = table_parser.add_argument_group("required arguments")
        table_group_required.add_argument("-i", "--input",           type=str,          required=True, nargs="*",     metavar="", help="Input file(s) and/or folder(s). '.tre' file(s) from ganon report.")
        table_group_required.add_argument("-e", "--input-extension", type=str,          default="tre", metavar="",                help="Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).")
        table_group_required.add_argument("-o", "--output-file",     type=str,          required=True,                            help="Output filename for the table")

        table_group_output = table_parser.add_argument_group("output arguments")
        table_group_output.add_argument("-l", "--output-value",  type=str,                    metavar="", default="counts", help="Output value on the table [percentage, counts]. percentage values are reported between [0-1]")
        table_group_output.add_argument("-f", "--output-format", type=str,                    metavar="", default="tsv",    help="Output format [tsv, csv]")
        table_group_output.add_argument("-t", "--top-sample",    type=unsigned_int(minval=0), metavar="", default=0,        help="Top hits of each sample individually")
        table_group_output.add_argument("-a", "--top-all",       type=unsigned_int(minval=0), metavar="", default=0,        help="Top hits of all samples (ranked by percentage)")
        table_group_output.add_argument("-m", "--min-frequency", type=int_or_float(minval=0), metavar="", default=0,        help="Minimum number/percentage of files containing an taxa to keep the taxa [values between 0-1 for percentage, >1 specific number]")
        table_group_output.add_argument("-r", "--rank",          type=str,   metavar="", default=None,   help="Define specific rank to report. Empty will report all ranks.")
        table_group_output.add_argument("-n", "--no-root",       action="store_true",    default=False,  help="Do not report root node entry and lineage. Direct and shared matches to root will be accounted as unclassified")
        table_group_output.add_argument("--header",              type=str,   metavar="", default="name", help="Header information [name, taxid, lineage]")
        table_group_output.add_argument("--unclassified-label",  type=str,   metavar="", default=None,   help="Add column with unclassified count/percentage with the chosen label. May be the same as --filtered-label (e.g. unassigned)")
        table_group_output.add_argument("--filtered-label",      type=str,   metavar="", default=None,   help="Add column with filtered count/percentage with the chosen label. May be the same as --unclassified-label (e.g. unassigned)")
        table_group_output.add_argument("--skip-zeros",          action="store_true",    default=False,  help="Do not print lines with only zero count/percentage")
        table_group_output.add_argument("--transpose",           action="store_true",    default=False,  help="Transpose output table (taxa as cols and files as rows)")

        table_group_optional = table_parser.add_argument_group("optional arguments")
        table_group_optional.add_argument("--verbose", action="store_true", default=False, help="Verbose output mode")
        table_group_optional.add_argument("--quiet",   action="store_true", default=False, help="Quiet output mode")

        ####################################################################################################

        filter_parser = argparse.ArgumentParser(add_help=False)
        filter_arguments = filter_parser.add_argument_group("filter arguments")
        filter_arguments.add_argument("--min-count",      type=int_or_float(minval=0), metavar="", default=0,  help="Minimum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1 specific number]")
        filter_arguments.add_argument("--max-count",      type=int_or_float(minval=0), metavar="", default=0,  help="Maximum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1 specific number]")
        filter_arguments.add_argument("--names",          type=str, nargs="*",         metavar="", default=[], help="Show only entries matching exact names of the provided list")
        filter_arguments.add_argument("--names-with",     type=str, nargs="*",         metavar="", default=[], help="Show entries containing full or partial names of the provided list")
        filter_arguments.add_argument("--taxids",         type=str, nargs="*",         metavar="", default=[], help="One or more taxids to report (including children taxa)")

        formatter_class = lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=120)
        subparsers = parser.add_subparsers()

        build = subparsers.add_parser("build",
                                      help="Download and build ganon default databases (refseq/genbank)",
                                      parents=[build_parser, build_default_parser, build_update_parser],
                                      formatter_class=formatter_class)
        build.set_defaults(which="build")

        build_custom = subparsers.add_parser("build-custom",
                                             help="Build custom ganon databases",
                                             parents=[build_custom_parser, build_default_parser, build_update_parser],
                                             formatter_class=formatter_class)
        build_custom.set_defaults(which="build-custom")

        update = subparsers.add_parser("update",
                                       help="Update ganon default databases",
                                       parents=[update_parser, build_update_parser],
                                       formatter_class=formatter_class)
        update.set_defaults(which="update")

        classify = subparsers.add_parser("classify",
                                         help="Classify reads against built databases",
                                         parents=[classify_parser],
                                         formatter_class=formatter_class)
        classify.set_defaults(which="classify")

        reassign = subparsers.add_parser("reassign",
                                         help="Reassign reads with multiple matches with an EM algorithm",
                                         parents=[reassign_parser],
                                         formatter_class=formatter_class)
        reassign.set_defaults(which="reassign")

        report = subparsers.add_parser("report",
                                       help="Generate reports from classification results",
                                       parents=[report_parser, filter_parser],
                                       formatter_class=formatter_class)
        report.set_defaults(which="report")

        table = subparsers.add_parser("table",
                                      help="Generate table from reports",
                                      parents=[table_parser, filter_parser],
                                      formatter_class=formatter_class)
        table.set_defaults(which="table")

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
        args = ["{}={}".format(k, repr(v)) for (k, v) in vars(self).items()]
        return "Config({})".format(", ".join(args))

    def set_defaults(self):
        if self.which in ["build", "build-custom"]:
            # If max-fp is not set, use default for ibf and hibf
            if self.max_fp is None:
                self.max_fp = 0.001 if self.filter_type == "hibf" else 0.05

        if self.which == "classify":
            if self.binning:
                self.rel_cutoff = [0.25]
                self.rel_filter = [0]
                self.min_count = 0
                self.report_type = "reads"
                
    def validate(self):

        if self.empty is True:
            print_log("Please provide one or more arguments")
            return False

        if self.verbose is True:
            self.quiet = False
        elif self.quiet is True:
            self.verbose = False

        if self.which == "build":

            if not self.organism_group and not self.taxid:
                print_log("--organism-group or --taxid required")
                return False

            if self.organism_group and self.taxid:
                print_log("--organism-group is mutually exclusive with --taxid")
                return False

        elif self.which == "build-custom":

            if self.input_file and self.input:
                print_log("--input-file is mutually exclusive with --input")
                return False

            if self.filter_type == "hibf" and self.hash_functions == 0:
                print_log("--filter-type hibf requires --hash-function value between 1 and 5")
                return False

            if self.level == "custom" and not self.input_file:
                print_log("--level custom requires --input-file")
                return False

            if self.level and self.level not in self.choices_level and self.taxonomy == "none":
                print_log("--taxonomy is required for --level " + self.level)
                return False

            if self.taxonomy == "ncbi":
                for entry in self.ncbi_sequence_info:
                    if entry not in self.choices_ncbi_sequence_info and not check_file(entry):
                        print_log("Invalid --get-sequence-info. Should be a valid file or: " +
                                  " ".join(self.choices_ncbi_sequence_info))
                        return False

                for entry in self.ncbi_file_info:
                    if entry not in self.choices_ncbi_file_info and not check_file(entry):
                        print_log("Invalid --get-file-info. Should be a valid file or: " +
                                  " ".join(self.choices_ncbi_file_info))
                        return False

        elif self.which == "update":
            if not check_folder(set_output_folder(self.db_prefix)):
                print_log("Folder to update not found: " + set_output_folder(self.db_prefix))
                return False

            if self.output_db_prefix and check_folder(set_output_folder(self.output_db_prefix)):
                print_log("New output folder already exists and it's not empty: " +
                          set_output_folder(self.output_db_prefix))
                return False

            update_required_files = ["_files/assembly_summary.txt",
                                     "_files/history.tsv",
                                     "_files/config.pkl"]
            for f in update_required_files:
                if not check_file(self.db_prefix + f):
                    print_log("Required file to update not found: " + f)
                    return False

        elif self.which == "classify":
            ibf = False
            hibf = False
            tax = 0
            for db_prefix in self.db_prefix:
                if check_file(db_prefix + ".hibf"):
                    hibf = True
                elif check_file(db_prefix + ".ibf"):
                    ibf = True
                else:
                    print_log("File not found: " + db_prefix + ".ibf/.hibf" )
                    return False

                if check_file(db_prefix + ".tax"):
                    tax += 1

            # Define use of HIBF and set hidden var
            if hibf and ibf:
                print_log(".ibf and .hibf filters cannot be used together in the same run" )
                return False
            elif hibf:
                # Hidden param
                self.hibf = True
                
            if tax < len(self.db_prefix) and tax > 0:
                print_log(".tax file has to be present for every .ibf/.hibf or none of them" )
                return False

            if not self.single_reads and not self.paired_reads:
                print_log("Please provide file[s] with --single-reads or --paired-reads")
                return False

            if self.single_reads:
                for f in self.single_reads:
                    if not check_file(f):
                        print_log("File not found: " + f )
                        return False

            if self.paired_reads:
                for f in self.paired_reads:
                    if not check_file(f):
                        print_log("File not found: " + f )
                        return False

                if len(self.paired_reads) % 2 != 0:
                    print_log("Invalid number of paired reads")
                    return False

            if not self.output_prefix and (self.output_all or self.output_one or self.output_unclassified):
                    print_log("--output-all / --output-one / --output-unclassified requires --output-prefix to be set")
                    return False

            if self.output_one and self.multiple_matches == "skip":
                print_log("--output-one requires --multiple-matches em/lca")
                return False              
        
        elif self.which == "report":

            if self.skip_hierarchy and self.keep_hierarchy:
                print_log("--skip-hierarchy and --keep-hierarchy are mutually exclusive")
                return False

            if self.db_prefix:
                for prefix in self.db_prefix:
                    file = prefix + ".tax" if not prefix.endswith(".tax") else prefix
                    if not check_file(file):
                        print_log("File not found: " + file)
                        return False

            if self.db_prefix and self.taxonomy == "skip":
                print_log("To skip taxonomy, omit --db-prefix and set --taxonomy skip")
                return False

        elif self.which == "table":
            pass

        return True

    def set_paths(self):
        missing_path = False
        self.ganon_path = self.ganon_path + "/" if self.ganon_path else ""

        if self.which in ["build", "update", "build-custom"]:
            ganon_genome_updater_paths = [self.ganon_path, self.ganon_path+"libs/genome_updater/", self.ganon_path+"../libs/genome_updater/"] if self.ganon_path else [None, "libs/genome_updater/"]
            for p in ganon_genome_updater_paths:
                self.path_exec["genome_updater"] = shutil.which("genome_updater.sh", path=p)
                if self.path_exec["genome_updater"] is not None: break
            if self.path_exec["genome_updater"] is None:
                print_log("genome_updater.sh was not found. Please inform a specific path with --ganon-path")
                missing_path = True

            ganon_get_seq_info_paths = [self.ganon_path, self.ganon_path+"scripts/", self.ganon_path+"../scripts/"] if self.ganon_path else [None, "scripts/"]
            for p in ganon_get_seq_info_paths:
                self.path_exec["get_seq_info"] = shutil.which("ganon-get-seq-info.sh", path=p)
                if self.path_exec["get_seq_info"] is not None: break
            if self.path_exec["get_seq_info"] is None:
                print_log("ganon-get-seq-info.sh script was not found. Please inform a specific path with --ganon-path")
                missing_path = True

            # if path is given, look for binaries only there
            ganon_build_paths = [self.ganon_path, self.ganon_path+"build/"] if self.ganon_path else [None, "build/"]
            for p in ganon_build_paths:
                self.path_exec["build"] = shutil.which("ganon-build", path=p)
                if self.path_exec["build"] is not None: break
            if self.path_exec["build"] is None:
                print_log("ganon-build binary was not found. Please inform a specific path with --ganon-path")
                missing_path = True

            if hasattr(self, 'filter_type') and self.filter_type == "hibf":
                self.raptor_path = self.raptor_path + "/" if self.raptor_path else ""
                raptor_paths = [self.raptor_path, self.raptor_path+"build/bin/"] if self.raptor_path else [None, "build/"]
                for p in raptor_paths:
                    self.path_exec["raptor"] = shutil.which("raptor", path=p)
                    if self.path_exec["raptor"] is not None: break
                if self.path_exec["raptor"] is None:
                    print_log("raptor binary was not found. Please inform a specific path with --raptor-path or use --filter-type ibf")
                    missing_path = True

        if self.which in ["classify"]:
            ganon_classify_paths = [self.ganon_path, self.ganon_path+"build/"] if self.ganon_path else [None, "build/"]
            for p in ganon_classify_paths:
                self.path_exec["classify"] = shutil.which("ganon-classify", path=p)
                if self.path_exec["classify"] is not None: break
            if self.path_exec["classify"] is None:
                print_log("ganon-classify binary was not found. Please inform a specific path with --ganon-path")
                missing_path = True

        return True if not missing_path else False


def unsigned_int(minval: int=0, maxval=None):
    def checker(val):
        val = int(val)
        if val < minval:
            raise argparse.ArgumentTypeError('%r must be >= %r' % (val, minval))
        if maxval is not None and val > maxval:
            raise argparse.ArgumentTypeError('%r must be <= %r' % (val, maxval))
        return val
    return checker


def unsigned_float(minval: int=0):
    def checker(val):
        val = float(val)
        if val is not None and val < minval:
            raise argparse.ArgumentTypeError('%r must be >= %r' % (val, minval))
        return val
    return checker


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
