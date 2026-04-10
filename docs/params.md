# Parameters

```
usage: ganon [-h] [-v]
             {build,build-custom,update,classify,reassign,report,table} ...

- - - - - - - - - -
   _  _  _  _  _   
  (_|(_|| |(_)| |  
   _|   v. 2.4.0
- - - - - - - - - -

positional arguments:
  {build,build-custom,update,classify,reassign,report,table}
    build               Download and build ganon default databases
                        (refseq/genbank)
    build-custom        Build custom ganon databases
    update              Update ganon default databases
    classify            Classify reads against built databases
    reassign            Reassign reads with multiple matches with an EM
                        algorithm
    report              Generate reports from classification results
    table               Generate table from reports

options:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
```

<details>
  <summary>ganon build</summary>

```
usage: ganon build [-h] [-g [ ...]] [-a [ ...]] [-l ] [-x ] [-m [ ...]] [-b [ ...]] [-o ] [-c] [-r] [-u ] [-z [ ...]]
                   [--skip-genome-size] [--download-threads ] -d DB_PREFIX [-t ] [-p ] [-k ] [-w ] [-s ] [-f ] [-j ]
                   [-y ] [-v ] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -g, --organism-group [ ...]
                        One or more organism groups to download [archaea, bacteria, fungi, human, invertebrate,
                        metagenomes, other, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral]. Mutually
                        exclusive --taxid (default: None)
  -a, --taxid [ ...]    One or more taxonomic identifiers to download. e.g. 562 (-x ncbi) or 's__Escherichia coli' (-x
                        gtdb). Mutually exclusive --organism-group (default: None)
  -d, --db-prefix DB_PREFIX
                        Database output prefix

database arguments:
  -l, --level           Highest level to build the database. Options: any available taxonomic rank [species, genus,
                        ...], 'leaves' for taxonomic leaves or 'assembly' for a assembly/strain based analysis (default:
                        species)

taxonomy arguments:
  -x, --taxonomy        Use taxonomy to enable taxonomic classification, lca and tax. reports [ncbi, gtdb, skip]
                        (default: ncbi)
  -m, --taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb: *taxonomy.tsv.gz (default: None)

download arguments:
  -b, --source [ ...]   Source to download [refseq, genbank] (default: ['refseq'])
  -o, --top             Download limited assemblies for each taxa. 0 for all. (default: 0)
  -c, --complete-genomes
                        Download only sub-set of complete genomes (default: False)
  -r, --reference-genomes
                        Download only sub-set of reference genomes (default: False)
  -u, --genome-updater 
                        Additional genome_updater parameters (https://github.com/pirovc/genome_updater) (default: None)
  -z, --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)
  --download-threads    Number of parallel sequence downloads from NCBI. (default: 8)

general arguments:
  -t, --threads         Number of sub-processes/threads to use (default: 1)
  -p, --max-fp          Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k, --kmer-size       The k-mer size to split sequences. (default: 19)
  -w, --window-size     The window-size to build filter with minimizers. (default: 31)
  -s, --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f, --filter-size     Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j, --mode            Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. Only valid for --filter-type ibf. (default:
                        avg)
  -y, --min-length      Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v, --filter-type     Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
                        path set with --raptor-path. --mode, --filter-size and --min-length will be ignored with hibf.
                        hibf will set --max-fp 0.001 as default. (default: hibf)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon build-custom</summary>

```
usage: ganon build-custom [-h] [-i [ ...]] [-e ] [-c] [-n ] [-a ] [-l ] [-z [ ...]] [--skip-genome-size] [-x ] [-b ]
                          [-m [ ...]] [-u [ ...]] [-g [ ...]] [--keep-invalid-taxa] [-r [ ...]] [-q [ ...]] -d DB_PREFIX
                          [-t ] [-p ] [-k ] [-w ] [-s ] [-f ] [-j ] [-y ] [-v ] [--restart] [--verbose] [--quiet]
                          [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -i, --input [ ...]    Input file(s) and/or folder(s). Mutually exclusive --input-file. (default: None)
  -e, --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: fna.gz)
  -c, --input-recursive
                        Look for files recursively in folder(s) provided with --input (default: False)
  -d, --db-prefix DB_PREFIX
                        Database output prefix

custom arguments:
  -n, --input-file      Tab-separated file with all necessary file/sequence information. Fields: file [<tab> target
                        <tab> node <tab> specialization <tab> specialization name]. For details:
                        https://pirovc.github.io/ganon/custom_databases/. Mutually exclusive --input (default: None)
  -a, --input-target    Target to use [file, sequence]. Parse input by file or by sequence. Using 'file' is recommended
                        and will speed-up the building process (default: file)
  -l, --level           Max. level to build the database. By default, --level is the --input-target. Options: any
                        available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy). Further
                        specialization options [assembly, custom]. assembly will retrieve and use the assembly accession
                        and name. custom requires and uses the specialization field in the --input-file. (default: None)
  -z, --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)

taxonomy arguments:
  -x, --taxonomy        Taxonomy matching the --input/--input-file. Enables taxonomic classification, lca and tax.
                        reports [ncbi, gtdb, gtdb-80, gtdb-83, gtdb-86.2, gtdb-89, gtdb-95, gtdb-202, gtdb-207,
                        gtdb-214.1, gtdb-220, gtdb-226, skip] (default: ncbi)
  -b, --convert-taxonomy 
                        Convert input taxonomy nodes (--taxonomy) to [ncbi-latest, gtdb-80, gtdb-83, gtdb-86.2, gtdb-89,
                        gtdb-95, gtdb-202, gtdb-207, gtdb-214.1, gtdb-220, gtdb-226]. (default: None)
  -m, --taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb: *taxonomy.tsv.gz (default: None)
  -u, --convert-taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi-latest: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb-version: *taxonomy.tsv.gz (default: None)
  -g, --convert-gtdb-files [ ...]
                        Use local gtdb conversion files instead of downloading. One for each version used in --taxonomy
                        and --convert-taxonomy. Files from https://github.com/pirovc/multitax/tree/main/data/gtdb
                        (default: None)
  --keep-invalid-taxa   Keep invalid taxa in the database, will be assigned to the root of the taxonomic tree. (default:
                        False)

ncbi arguments:
  -r, --ncbi-sequence-info [ ...]
                        Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information.
                        [eutils, nucl_gb, nucl_wgs, nucl_est, nucl_gss, pdb, prot, dead_nucl, dead_wgs, dead_prot or one
                        or more accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/].
                        By default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise. (default:
                        [])
  -q, --ncbi-file-info [ ...]
                        Downloads assembly_summary files to extract target information. [refseq, genbank,
                        refseq_historical, genbank_historical or one or more assembly_summary files from
                        https://ftp.ncbi.nlm.nih.gov/genomes/] (default: ['refseq', 'genbank'])

general arguments:
  -t, --threads         Number of sub-processes/threads to use (default: 1)
  -p, --max-fp          Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k, --kmer-size       The k-mer size to split sequences. (default: 19)
  -w, --window-size     The window-size to build filter with minimizers. (default: 31)
  -s, --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f, --filter-size     Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j, --mode            Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. Only valid for --filter-type ibf. (default:
                        avg)
  -y, --min-length      Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v, --filter-type     Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
                        path set with --raptor-path. --mode, --filter-size and --min-length will be ignored with hibf.
                        hibf will set --max-fp 0.001 as default. (default: hibf)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon update</summary>

```
usage: ganon update [-h] -d DB_PREFIX [-o ] [-t ] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -d, --db-prefix DB_PREFIX
                        Existing database input prefix

general arguments:
  -o, --output-db-prefix 
                        Output database prefix. By default will be the same as --db-prefix and overwrite files (default:
                        None)
  -t, --threads         Number of sub-processes/threads to use (default: 1)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon classify</summary>

```
usage: ganon classify [-h] -d [DB_PREFIX ...] -o OUTPUT_PREFIX [-s [reads.fq[.gz] ...]]
                      [-p [reads.1.fq[.gz] reads.2.fq[.gz] ...]] [-a [file.tsv ...]] [-c [ ...]] [-e [ ...]] [-m ]
                      [--ranks [ ...]] [--min-count ] [--report-type ] [--skip-report] [--output-one] [--output-all]
                      [--output-unclassified] [--output-stats] [--output-single] [-t ] [-b] [-f [ ...]] [-l [ ...]]
                      [--verbose] [--quiet]

options:
  -h, --help            show this help message and exit

required arguments:
  -d, --db-prefix [DB_PREFIX ...]
                        Database input prefix[es]
  -o, --output-prefix OUTPUT_PREFIX
                        Output prefix for base report (.rep) and tree-like report (.tre).
  -s, --single-reads [reads.fq[.gz] ...]
                        Multi-fastq[.gz] file[s] to classify (default: None)
  -p, --paired-reads [reads.1.fq[.gz] reads.2.fq[.gz] ...]
                        Multi-fastq[.gz] pairs of file[s] to classify (default: None)
  -a, --batch-reads [file.tsv ...]
                        File with single- or paired-end reads to be processed in one run. Prefix can be repeated.
                        Example: prefix <tab> file1 [<tab> file2] (default: None)

cutoff/filter arguments:
  -c, --rel-cutoff [ ...]
                        Min. percentage of a read (set of k-mers) shared with a reference necessary to consider a match.
                        Generally used to remove low similarity matches. Single value or one per database (e.g. 0.7 1
                        0.25). 0 for no cutoff (default: [0.75])
  -e, --rel-filter [ ...]
                        Additional relative percentage of matches (relative to the best match) to keep. Generally used
                        to keep top matches above cutoff. Single value or one per hierarchy (e.g. 0.1 0). 1 for no
                        filter (default: [0.1])

post-processing/report arguments:
  -m, --multiple-matches 
                        Method to solve reads with multiple matches [em, lca, skip]. em -> expectation maximization
                        algorithm based on unique matches. lca -> lowest common ancestor based on taxonomy. The EM
                        algorithm can be executed later with 'ganon reassign' using the .all file (--output-all).
                        (default: em)
  --ranks [ ...]        Ranks to report taxonomic abundances (.tre). empty will report default ranks [domain phylum
                        class order family genus species assembly]. (default: [])
  --min-count           Minimum percentage/counts to report an taxa (.tre) [use values between 0-1 for percentage, >1
                        for counts] (default: 5e-05)
  --report-type         Type of report (.tre) [abundance, reads, matches, dist, corr]. More info in 'ganon report'.
                        (default: abundance)
  --skip-report         Disable tree-like report (.tre) at the end of classification. Can be done later with 'ganon
                        report'. (default: False)

output arguments:
  --output-one          Output a file with one match for each read (.one) either an unique match or a result from the EM
                        or a LCA algorithm (--multiple-matches) (default: False)
  --output-all          Output a file with all unique and multiple matches (.all) (default: False)
  --output-unclassified
                        Output a file with unclassified read headers (.unc) (default: False)
  --output-stats        Output a file with statistic of classification (.sta) (default: False)
  --output-single       When using multiple hierarchical levels, output everything in one file instead of one per
                        hierarchy (default: False)

other arguments:
  -t, --threads         Number of sub-processes/threads to use (default: 1)
  -b, --binning         Optimized parameters for binning (--rel-cutoff 0.25 --rel-filter 0 --min-count 0 --report-type
                        reads). Will report sequence abundances (.tre) instead of tax. abundance. (default: False)
  -f, --fpr-query [ ...]
                        Max. false positive of a query to accept a match. Applied after --rel-cutoff and --rel-filter.
                        Generally used to remove false positives matches querying a database build with large --max-fp.
                        Single value or one per hierarchy (e.g. 0.1 0). 1 for no filter (default: [1e-05])
  -l, --hierarchy-labels [ ...]
                        Hierarchy definition of --db-prefix files to be classified. Can also be a string, but input will
                        be sorted to define order (e.g. 1 1 2 3). The default value reported without hierarchy is 'H1'
                        (default: None)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
```

</details>

<details>
  <summary>ganon reassign</summary>

```
usage: ganon reassign [-h] -i [ ...] [-o OUTPUT_PREFIX] [-e ] [-s ] [--remove-all] [--skip-one] [--skip-rep] [--verbose]
                      [--quiet]

options:
  -h, --help            show this help message and exit

required arguments:
  -i, --input-prefix [ ...]
                        Input prefix to find files from ganon classify (.rep and .all)
  -o, --output-prefix OUTPUT_PREFIX
                        Alternative output prefix for reassigned files. If not provided, will use same path of input
                        files (will overwrite .rep). In case of multiple files, the output will be the suffix. Example:
                        {output_prefix}{filename}.one (default: )

EM arguments:
  -e, --max-iter        Max. number of iterations for the EM algorithm. If 0, will run until convergence (check
                        --threshold) (default: 10)
  -s, --threshold       Convergence threshold limit to stop the EM algorithm. (default: 0)

other arguments:
  --remove-all          Remove input file (.all) after processing. (default: False)
  --skip-one            Do not write output file (.one) after processing. (default: False)
  --skip-rep            Do not write report file (.rep) after processing. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
```

</details>

<details>
  <summary>ganon report</summary>

```
usage: ganon report [-h] -i [ ...] [-e INPUT_EXTENSION] [-d [ ...]] [-x ] [-m [ ...]] [-z [ ...]] [--skip-genome-size]
                    [-o OUTPUT_PREFIX] [-f ] [-t ] [-r [ ...]] [-s ] [-a] [-y] [-p [ ...]] [-k [ ...]] [-c ] [-n]
                    [--verbose] [--quiet] [--min-count ] [--max-count ] [--names [ ...]] [--names-with [ ...]]
                    [--taxids [ ...]]

options:
  -h, --help            show this help message and exit

required arguments:
  -i, --input [ ...]    Input file(s) and/or folder(s). '.rep' file(s) from ganon classify.
  -e, --input-extension INPUT_EXTENSION
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: rep)

db/tax arguments:
  -d, --db-prefix [ ...]
                        Database prefix(es) used for classification. Only '.tax' file(s) are required. If not provided,
                        new taxonomy will be downloaded. Mutually exclusive with --taxonomy. (default: [])
  -x, --taxonomy        Taxonomy database to use [ncbi, gtdb, skip]. Mutually exclusive with --db-prefix. (default:
                        ncbi)
  -m, --taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb: *taxonomy.tsv.gz (default: None)
  -z, --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Valid only without --db-prefix. Activate this option when
                        using sequences not representing full genomes. (default: False)

output arguments:
  -o, --output-prefix OUTPUT_PREFIX
                        Output prefix for report file 'output_prefix.tre'. In case of multiple files, the base input
                        filename will be appended at the end of the output file 'output_prefix + FILENAME.tre' (default:
                        )
  -f, --output-format   Output format [text, tsv, csv, bioboxes]. text outputs a tabulated formatted text file for
                        better visualization. bioboxes is the the CAMI challenge profiling format (only
                        percentage/abundances are reported). (default: tsv)
  -t, --report-type     Type of report [abundance, reads, matches, dist, corr]. 'abundance' -> tax. abundance (re-
                        distribute read counts and correct by genome size), 'reads' -> sequence abundance, 'matches' ->
                        report all unique and shared matches, 'dist' -> like reads with re-distribution of shared read
                        counts only, 'corr' -> like abundance without re-distribution of shared read counts (default:
                        abundance)
  -r, --ranks [ ...]    Ranks to report ['', 'all', custom list]. 'all' for all possible ranks. empty for default ranks
                        [domain phylum class order family genus species assembly]. (default: [])
  -s, --sort            Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage
                        (with --ranks all) (default: )
  -a, --no-orphan       Omit orphan nodes from the final report. Otherwise, orphan nodes (= nodes not found in the
                        db/tax) are reported as 'na' with root as direct parent. (default: False)
  -y, --split-hierarchy
                        Split output reports by hierarchy (from ganon classify --hierarchy-labels). If activated, the
                        output files will be named as '{output_prefix}.{hierarchy}.tre' (default: False)
  -p, --skip-hierarchy [ ...]
                        One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels) (default:
                        [])
  -k, --keep-hierarchy [ ...]
                        One or more hierarchies to keep in the report (from ganon classify --hierarchy-labels) (default:
                        [])
  -c, --top-percentile 
                        Top percentile filter, based on percentage/relative abundance. Applied only at default ranks
                        [domain phylum class order family genus species assembly] (default: 0)
  -n, --normalize       Ignore the number of unclassified reads, normalizing the output to 100%. Use with caution, can
                        drastically change abundance estimations. (default: False)

optional arguments:
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)

filter arguments:
  --min-count           Minimum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --max-count           Maximum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --names [ ...]        Show only entries matching exact names of the provided list (default: [])
  --names-with [ ...]   Show entries containing full or partial names of the provided list (default: [])
  --taxids [ ...]       One or more taxids to report (including children taxa) (default: [])
```

</details>

<details>
  <summary>ganon table</summary>

```
usage: ganon table [-h] -i [ ...] [-e ] -o OUTPUT_FILE [-l ] [-f ] [-t ] [-a ] [-m ] [-r ] [-n] [--header ]
                   [--unclassified-label ] [--filtered-label ] [--skip-zeros] [--transpose] [--verbose] [--quiet]
                   [--min-count ] [--max-count ] [--names [ ...]] [--names-with [ ...]] [--taxids [ ...]]

options:
  -h, --help            show this help message and exit

required arguments:
  -i, --input [ ...]    Input file(s) and/or folder(s). '.tre' file(s) from ganon report.
  -e, --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: tre)
  -o, --output-file OUTPUT_FILE
                        Output filename for the table

output arguments:
  -l, --output-value    Output value on the table [percentage, counts]. percentage values are reported between [0-1]
                        (default: counts)
  -f, --output-format   Output format [tsv, csv] (default: tsv)
  -t, --top-sample      Top hits of each sample individually (default: 0)
  -a, --top-all         Top hits of all samples (ranked by percentage) (default: 0)
  -m, --min-frequency   Minimum number/percentage of files containing an taxa to keep the taxa [values between 0-1 for
                        percentage, >1 specific number] (default: 0)
  -r, --rank            Define specific rank to report. Empty will report all ranks. (default: None)
  -n, --no-root         Do not report root node entry and lineage. Direct and shared matches to root will be accounted
                        as unclassified (default: False)
  --header              Header information [name, taxid, lineage] (default: name)
  --unclassified-label 
                        Add column with unclassified count/percentage with the chosen label. May be the same as
                        --filtered-label (e.g. unassigned) (default: None)
  --filtered-label      Add column with filtered count/percentage with the chosen label. May be the same as
                        --unclassified-label (e.g. unassigned) (default: None)
  --skip-zeros          Do not print lines with only zero count/percentage (default: False)
  --transpose           Transpose output table (taxa as cols and files as rows) (default: False)

optional arguments:
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)

filter arguments:
  --min-count           Minimum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --max-count           Maximum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --names [ ...]        Show only entries matching exact names of the provided list (default: [])
  --names-with [ ...]   Show entries containing full or partial names of the provided list (default: [])
  --taxids [ ...]       One or more taxids to report (including children taxa) (default: [])
```

</details>
