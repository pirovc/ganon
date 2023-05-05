# ganon

[![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1101%2F406017-blue)](https://dx.doi.org/10.1093/bioinformatics/btaa458)

ganon is designed to index large sets of genomic reference sequences and to classify reads against them efficiently. The tool uses Interleaved Bloom Filters as indices based on k-mers/minimizers. It was mainly developed, but not limited, to the metagenomics classification problem: quickly assign sequence fragments to their closest reference among thousands of references. After classification, taxonomic abundance is estimated and reported.

## Features

- NCBI and GTDB native support for taxonomic classification (also runs without taxonomy)
- integrated download of commonly used reference sequences from RefSeq/Genbank (`ganon build`)
- update indices incrementally (`ganon update`)
- customizable build for pre-downloaded or non-standard sequence files (`ganon build-custom`)
- build and classify at different taxonomic levels, file, sequence, strain/assembly or custom specialization
- perform [hierarchical classification](#multiple-and-hierarchical-classification): use several databases in any order
- [report](#report) the lowest common ancestor (LCA) but also multiple and unique matches for every read
- [report](#report) sequence or taxonomic abundances as well as total number of matches
- reassignment of reads with multiple matches to a unique match with an EM algorithm
- generate reports and contingency tables for multi-sample studies with several filter options

ganon achieved very good results in [our own evaluations](https://dx.doi.org/10.1093/bioinformatics/btaa458) but also in independent evaluations: [LEMMI](https://lemmi-v1.ezlab.org/), [LEMMI v2](https://lemmi.ezlab.org/) and [CAMI2](https://dx.doi.org/10.1038/s41592-022-01431-4)

## Installation with conda

The easiest way to install ganon is via conda, using the bioconda and conda-forge channels:

```bash
conda install -c bioconda -c conda-forge ganon
```

However, there are possible performance benefits compiling ganon from source in the target machine rather than using the conda version. To do so, please follow the instructions below:

## Installation from source

### Dependencies

- gcc >=7
- cmake >=3.10
- zlib
- python >=3.6
- pandas >=1.1.0
- [multitax](https://github.com/pirovc/multitax) >=1.3.1

```bash
python3 -V # >=3.6
python3 -m pip install "pandas>=1.1.0" "multitax>=1.3.1"
```

### Downloading and building ganon + submodules

```bash
git clone --recurse-submodules https://github.com/pirovc/ganon.git
```
  
```bash
cd ganon
python3 setup.py install --record files.txt  # optional
mkdir build_cpp
cd build_cpp
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF -DLONGREADS=OFF ..
make
sudo make install  # optional
```

- to change install location (e.g. `/myprefix/bin/`), set the installation prefix in the cmake command with `-DCMAKE_INSTALL_PREFIX=/myprefix/ `
- use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.
- to classify extremely large reads (>65000bp) use `-DLONGREADS=ÒN`

If everything was properly installed, the following commands should show the help pages without errors:

```bash
ganon -h
```

### Running tests

```bash
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/  # optional - downloads large files
cd build_cpp/
ctest -VV .
```
 
## Important parameters

The most important parameters and trade-offs are:

- `ganon build` `--hibf`: build smaller databases that can be queried faster. Building will take longer.
- `ganon build` `--window-size --kmer-size`: the *window* value should always be the same or larger than the *kmer* value. The larger the difference between them, the smaller the database will be. However, some sensitivity/precision loss in classification is expected with small *kmer* and/or large *window*. Larger *kmer* values (e.g. `31`) will improve classification, specially read binning, at a cost of way bigger databases.
---
- `ganon classify` `--rel-cutoff`: this value defines the threshold for matches between reads and database. Higher `--rel-cutoff` values will improve precision and decrease sensitivity with expected less unique matches but an increase in overall matches. For taxonomic profiling, a higher value between `0.4` and `0.8` may provide better results. For read binning, lower values between `0.2` and `0.4` are recommended.
- `ganon classify` `--rel-filter`: further filter top matches after cutoff is applied. Usually set between `0` and `0.2`.
- `ganon classify` `--reassign`: runs an EM-algorithm to reassign reads that received multiple matches. It provides a unique match for each read at the level the database was built (e.g. assembly or species). Mostly useful for read binning, with little overall impact on taxonomic profiling. Can be used independently with `ganon reassign`.
---
- `ganon report` `--report-type`: reports either taxonomic, sequence or matches abundances. Use `corr` or `abundance` for taxonomic profiling, `reads` or `dist` for sequence profiling and `matches` to report a summary of all matches.
- `ganon report` `--min-count`: cutoff to discard underrepresented taxa. Useful to remove the common long tail of spurious matches and false positives when performing classification. Values between `0.0001` (0.01%) and `0.001` (0.1%) improved sensitivity and precision in our evaluations. The higher the value, the more precise the outcome, with a sensitivity loss. Alternatively `--top-percentile` can be used to keep a relative amount of taxa instead a hard cutoff.

The numeric values above are averages from several experiments with different sample types and database contents. They may not work as expected for your data. If you are not sure which values to use or see something unexpected, please open an [issue](https://github.com/pirovc/ganon/issues).

## Parameters

```
usage: ganon [-h] [-v]
             {build,build-custom,update,classify,reassign,report,table} ...

- - - - - - - - - -
   _  _  _  _  _   
  (_|(_|| |(_)| |  
   _|   v. 1.5.0
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
usage: ganon build [-h] [-g [...]] [-a [...]] [-b [...]] [-o] [-c] [-u] [-m [...]] [-z [...]] -d DB_PREFIX [-x] [-t]
                   [-p] [-f] [-k] [-w] [-s] [-j] [--hibf] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -g [ ...], --organism-group [ ...]
                        One or more organism groups to download [archaea, bacteria, fungi, human, invertebrate,
                        metagenomes, other, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral]. Mutually
                        exclusive --taxid (default: None)
  -a [ ...], --taxid [ ...]
                        One or more taxonomic identifiers to download. e.g. 562 (-x ncbi) or 's__Escherichia coli' (-x
                        gtdb). Mutually exclusive --organism-group (default: None)
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (default: None)

download arguments:
  -b [ ...], --source [ ...]
                        Source to download [refseq, genbank] (default: ['refseq'])
  -o , --top            Download limited assemblies for each taxa. 0 for all. (default: 0)
  -c, --complete-genomes
                        Download only sub-set of complete genomes (default: False)
  -u , --genome-updater 
                        Additional genome_updater parameters (https://github.com/pirovc/genome_updater) (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)

important arguments:
  -x , --taxonomy       Set taxonomy to enable taxonomic classification, lca and reports [ncbi, gtdb, skip] (default:
                        ncbi)
  -t , --threads 

advanced arguments:
  -p , --max-fp         Max. false positive rate for bloom filters Mutually exclusive --filter-size. (default: 0.05)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. (default: 0)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [0-5]. 0 to detect optimal value.
                        (default: 4)
  -j , --mode           Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. (default: avg)
  --hibf                Builds an HIBF with raptor/chopper (v3). --mode and --filter-size will be ignored. (default:
                        False)

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
usage: ganon build-custom [-h] [-i [...]] [-e] [-n] [-a] [-l] [-m [...]] [-z [...]] [-r [...]] [-q [...]] -d DB_PREFIX
                          [-x] [-t] [-p] [-f] [-k] [-w] [-s] [-j] [--hibf] [--restart] [--verbose] [--quiet]
                          [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). Mutually exclusive --input-file. (default: None)
  -e , --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: fna.gz)
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (default: None)

custom arguments:
  -n , --input-file     Manually set information for input files: file <tab> [target <tab> node <tab> specialization
                        <tab> specialization name]. target is the sequence identifier if --input-target sequence (file
                        can be repeated for multiple sequences). if --input-target file and target is not set, filename
                        is used. node is the taxonomic identifier. Mutually exclusive --input (default: None)
  -a , --input-target   Target to use [file, sequence]. By default: 'file' if multiple input files are provided or
                        --input-file is set, 'sequence' if a single file is provided. Using 'file' is recommended and
                        will speed-up the building process (default: None)
  -l , --level          Use a specialized target to build the database. By default, --level is the --input-target.
                        Options: any available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy).
                        Further specialization options [assembly, custom]. assembly will retrieve and use the assembly
                        accession and name. custom requires and uses the specialization field in the --input-file.
                        (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)

ncbi arguments:
  -r [ ...], --ncbi-sequence-info [ ...]
                        Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information.
                        [eutils, nucl_gb, nucl_wgs, nucl_est, nucl_gss, pdb, prot, dead_nucl, dead_wgs, dead_prot or one
                        or more accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/].
                        By default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise. (default:
                        [])
  -q [ ...], --ncbi-file-info [ ...]
                        Downloads assembly_summary files to extract target information. [refseq, genbank,
                        refseq_historical, genbank_historical or one or more assembly_summary files from
                        https://ftp.ncbi.nlm.nih.gov/genomes/] (default: ['refseq', 'genbank'])

important arguments:
  -x , --taxonomy       Set taxonomy to enable taxonomic classification, lca and reports [ncbi, gtdb, skip] (default:
                        ncbi)
  -t , --threads 

advanced arguments:
  -p , --max-fp         Max. false positive rate for bloom filters Mutually exclusive --filter-size. (default: 0.05)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. (default: 0)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [0-5]. 0 to detect optimal value.
                        (default: 4)
  -j , --mode           Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. (default: avg)
  --hibf                Builds an HIBF with raptor/chopper (v3). --mode and --filter-size will be ignored. (default:
                        False)

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
usage: ganon update [-h] -d DB_PREFIX [-o] [-t] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Existing database input prefix (default: None)

important arguments:
  -o , --output-db-prefix 
                        Output database prefix. By default will be the same as --db-prefix and overwrite files (default:
                        None)
  -t , --threads 

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
usage: ganon classify [-h] -d [DB_PREFIX ...] [-s [reads.fq[.gz] ...]] [-p [reads.1.fq[.gz] reads.2.fq[.gz] ...]]
                      [-c [...]] [-e [...]] [-o] [--output-lca] [--output-all] [--output-unclassified] [--output-single]
                      [-t] [-l [...]] [-r [...]] [-a] [--verbose] [--quiet]

options:
  -h, --help            show this help message and exit

required arguments:
  -d [DB_PREFIX ...], --db-prefix [DB_PREFIX ...]
                        Database input prefix[es] (default: None)
  -s [reads.fq[.gz] ...], --single-reads [reads.fq[.gz] ...]
                        Multi-fastq[.gz] file[s] to classify (default: None)
  -p [reads.1.fq[.gz] reads.2.fq[.gz] ...], --paired-reads [reads.1.fq[.gz] reads.2.fq[.gz] ...]
                        Multi-fastq[.gz] pairs of file[s] to classify (default: None)

cutoff/filter arguments:
  -c [ ...], --rel-cutoff [ ...]
                        Min. percentage of a read (set of minimizers) shared with the a reference necessary to consider
                        a match. Generally used to cutoff low similarity matches. Single value or one per database (e.g.
                        0.7 1 0.25). 0 for no cutoff (default: [0.75])
  -e [ ...], --rel-filter [ ...]
                        Additional relative percentage of minimizers (relative to the best match) to keep a match.
                        Generally used to select best matches above cutoff. Single value or one per hierarchy (e.g. 0.1
                        0). 1 for no filter (default: [0.0])

output arguments:
  -o , --output-prefix 
                        Output prefix for output (.rep) and report (.tre). Empty to output to STDOUT (only .rep)
                        (default: None)
  --output-lca          Output an additional file with one lca match for each read (.lca) (default: False)
  --output-all          Output an additional file with all matches. File can be very large (.all) (default: False)
  --output-unclassified
                        Output an additional file with unclassified read headers (.unc) (default: False)
  --output-single       When using multiple hierarchical levels, output everything in one file instead of one per
                        hierarchy (default: False)

other arguments:
  -t , --threads        Number of sub-processes/threads to use (default: 1)
  -l [ ...], --hierarchy-labels [ ...]
                        Hierarchy definition of --db-prefix files to be classified. Can also be a string, but input will
                        be sorted to define order (e.g. 1 1 2 3). The default value reported without hierarchy is 'H1'
                        (default: None)
  -r [ ...], --ranks [ ...]
                        Ranks to report taxonomic abundances (.tre). empty will report default ranks [superkingdom,
                        phylum, class, order, family, genus, species, assembly]. This file can be re-generated with the
                        'ganon report' command for other types of abundances (reads, matches) with further filtration
                        and output options (default: [])
  -a, --reassign        Reassign reads with multiple matches with an EM algorithm. Will enforce --output-all. This file
                        can be re-generated with the 'ganon reassign'. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
```

</details>

<details>
  <summary>ganon reassign</summary>

```
usage: ganon reassign [-h] -i  -o OUTPUT_PREFIX [-e] [-s] [--verbose] [--quiet]

options:
  -h, --help            show this help message and exit

required arguments:
  -i , --input-prefix   Input prefix to find files from ganon classify (.all and optionally .rep) (default: None)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix for reassigned file (.all and optionally .rep). In case of multiple files, the
                        base input filename will be appended at the end of the output file 'output_prefix +
                        FILENAME.all' (default: None)

EM arguments:
  -e , --max-iter       Max. number of iterations for the EM algorithm. If 0, will run until convergence (check
                        --threshold) (default: 10)
  -s , --threshold      Convergence threshold limit to stop the EM algorithm. (default: 0)

other arguments:
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
```

</details>

<details>
  <summary>ganon report</summary>

```
usage: ganon report [-h] -i [...] [-e INPUT_EXTENSION] -o OUTPUT_PREFIX [-d [...]] [-x] [-m [...]] [-z [...]] [-f] [-t]
                    [-r [...]] [-s] [-a] [-y] [-p [...]] [-k [...]] [-c] [--verbose] [--quiet] [--min-count]
                    [--max-count] [--names [...]] [--names-with [...]] [--taxids [...]]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). '.rep' file(s) from ganon classify. (default: None)
  -e INPUT_EXTENSION, --input-extension INPUT_EXTENSION
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: rep)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix for report file 'output_prefix.tre'. In case of multiple files, the base input
                        filename will be appended at the end of the output file 'output_prefix + FILENAME.tre' (default:
                        None)

db/tax arguments:
  -d [ ...], --db-prefix [ ...]
                        Database prefix(es) used for classification. Only '.tax' file(s) are required. If not provided,
                        new taxonomy will be downloaded. Mutually exclusive with --taxonomy. (default: [])
  -x , --taxonomy       Taxonomy database to use [ncbi, gtdb, skip]. Mutually exclusive with --db-prefix. (default:
                        ncbi)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)

output arguments:
  -f , --output-format 
                        Output format [text, tsv, csv, bioboxes]. text outputs a tabulated formatted text file for
                        better visualization. bioboxes is the the CAMI challenge profiling format (only
                        percentage/abundances are reported). (default: tsv)
  -t , --report-type    Type of report [abundance, reads, matches, dist, corr]. 'abundance' -> tax. abundance (re-
                        distribute read counts and correct by genome size), 'reads' -> sequence abundance, 'matches' ->
                        report all unique and shared matches, 'dist' -> like reads with re-distribution of shared read
                        counts only, 'corr' -> like abundance without re-distribution of shared read counts (default:
                        abundance)
  -r [ ...], --ranks [ ...]
                        Ranks to report ['', 'all', custom list]. 'all' for all possible ranks. empty for default ranks
                        [superkingdom, phylum, class, order, family, genus, species, assembly]. (default: [])
  -s , --sort           Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage
                        (with --ranks all) (default: )
  -a, --no-orphan       Omit orphan nodes from the final report. Otherwise, orphan nodes (= nodes not found in the
                        db/tax) are reported as 'na' with root as direct parent. (default: False)
  -y, --split-hierarchy
                        Split output reports by hierarchy (from ganon classify --hierarchy-labels). If activated, the
                        output files will be named as '{output_prefix}.{hierarchy}.tre' (default: False)
  -p [ ...], --skip-hierarchy [ ...]
                        One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels) (default:
                        [])
  -k [ ...], --keep-hierarchy [ ...]
                        One or more hierarchies to keep in the report (from ganon classify --hierarchy-labels) (default:
                        [])
  -c , --top-percentile 
                        Top percentile filter, based on percentage/relative abundance. Applied only at default ranks
                        [superkingdom, phylum, class, order, family, genus, species, assembly] (default: 0)

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
usage: ganon table [-h] -i [...] [-e] -o OUTPUT_FILE [-l] [-f] [-t] [-a] [-m] [-r] [-n] [--header]
                   [--unclassified-label] [--filtered-label] [--skip-zeros] [--transpose] [--verbose] [--quiet]
                   [--min-count] [--max-count] [--names [...]] [--names-with [...]] [--taxids [...]]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). '.tre' file(s) from ganon report. (default: None)
  -e , --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: tre)
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output filename for the table (default: None)

output arguments:
  -l , --output-value   Output value on the table [percentage, counts]. percentage values are reported between [0-1]
                        (default: counts)
  -f , --output-format 
                        Output format [tsv, csv] (default: tsv)
  -t , --top-sample     Top hits of each sample individually (default: 0)
  -a , --top-all        Top hits of all samples (ranked by percentage) (default: 0)
  -m , --min-frequency 
                        Minimum number/percentage of files containing an taxa to keep the taxa [values between 0-1 for
                        percentage, >1 specific number] (default: 0)
  -r , --rank           Define specific rank to report. Empty will report all ranks. (default: None)
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
