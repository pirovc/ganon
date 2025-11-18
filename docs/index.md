# ganon [![GitHub release (latest by date)](https://img.shields.io/github/v/release/pirovc/ganon)](https://github.com/pirovc/ganon)

[![Build Status](https://app.travis-ci.com/pirovc/ganon.svg?token=q6Nfx8pLHh8hV3hLz3Pq&branch=master)](https://app.travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1101%2F406017-blue)](https://dx.doi.org/10.1093/bioinformatics/btaa458) 

Code: [GitHub repository](https://github.com/pirovc/ganon)

[ganon2 pre-print](https://www.biorxiv.org/content/10.1101/2023.12.07.570547)

ganon is designed to index large sets of genomic reference sequences and to classify reads against them efficiently. The tool uses [Hierarchical Interleaved Bloom Filters](https://doi.org/10.1186/s13059-023-02971-4) as indices based on k-mers with optional minimizers. It was mainly developed, but not limited, to the metagenomics classification problem: quickly assign sequence fragments to their closest reference among thousands of references. After classification, taxonomic or sequence abundances are estimated and reported.

## Features

- integrated download and build of any subset from [RefSeq/Genbank/GTDB](default_databases.md#refseq-and-genbank) with incremental [updates](default_databases.md#update-ganon-update)
- NCBI and [GTDB](default_databases.md#gtdb) native support for taxonomic classification, custom taxonomy or no taxonomy at all
- [customizable database](custom_databases.md) build for local or non-standard sequence files
- optimized [taxonomic binning](classification.md#binning) and [profiling](classification.md#profiling) configurations
- build and classify at various taxonomic levels, strain, assembly, file, sequence or custom specialization
- [hierarchical classification](classification.md#multiple-and-hierarchical-classification) using several databases in one or more levels in just one run
- [EM and/or LCA](classification.md#reads-with-multiple-matches) algorithms to solve multiple-matching reads
- reporting of multiple and unique matches for every read
- [reporting](reports.md#report-type-report-type) of sequence, taxonomic or multi-match abundances with optional genome size correction
- advanced tree-like [reports](reports.md) with several filter options
- generation of [contingency tables](table.md) with several filters for multi-sample studies

ganon achieved very good results in [our own evaluations](https://dx.doi.org/10.1093/bioinformatics/btaa458) but also in independent evaluations: [LEMMI](https://lemmi-v1.ezlab.org/), [LEMMI v2](https://lemmi.ezlab.org/) and [CAMI2](https://dx.doi.org/10.1038/s41592-022-01431-4)

## Installation with conda

The easiest way to install ganon is via conda, using the bioconda and conda-forge channels:

```bash
conda install -c bioconda -c conda-forge ganon
```

However, there are possible performance benefits compiling ganon from source in the target machine rather than using the conda version. To do so, please follow the instructions below:

## Installation from source

### Python dependencies

- python >=3.6
- pandas >=1.2.0
- [multitax](https://github.com/pirovc/multitax) >=1.3.1
- [genome_updater](https://github.com/pirovc/genome_updater) >=0.6.4

```bash
# Python version should be >=3.6
python3 -V

# Install packages via pip or conda:
# PIP
python3 -m pip install "pandas>=1.2.0" "multitax>=1.3.1"
wget --quiet --show-progress https://raw.githubusercontent.com/pirovc/genome_updater/master/genome_updater.sh && chmod +x genome_updater.sh

# Conda/Mamba (alternative)
conda install -c bioconda -c conda-forge "pandas>=1.2.0" "multitax>=1.3.1" "genome_updater>=0.6.4"
```
### C++ dependencies

- GCC >=11
- CMake >=3.5
- zlib
- bzip2
- raptor ==3.0.1

!!! tip
    If your system has GCC version 10 or below, you can create an environment with the latest conda-forge GCC version and dependencies: `conda create -c conda-forge -n gcc-conda cxx-compiler zlib bzip2 "cmake>=3.5"` and activate the environment with: `source activate gcc-conda`.
    
    In CMake, you may have set the environment include directory with the following parameter: `-DSEQAN3_CXX_FLAGS="-I/path/to/miniconda3/envs/gcc-conda/include/"` changing `/path/to/miniconda3` with your local path to the conda installation.

### Downloading and building ganon + submodules

```bash
git clone --recurse-submodules https://github.com/pirovc/ganon.git
```
  
```bash
# Install Python side
cd ganon
pip install .

# Compile and install C++ side
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF -DLONGREADS=OFF ..
make -j 4
sudo make install  # optional
```

- to change install location (e.g. `/myprefix/bin/`), set the installation prefix in the cmake command with `-DCMAKE_INSTALL_PREFIX=/myprefix/ `
- use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.
- to classify extremely large reads or contigs that would need more than 65000 k-mers, use `-DLONGREADS=ON`

### Installing raptor

The easiest way to install [raptor](https://github.com/seqan/raptor) is via conda with `conda install -c bioconda -c conda-forge "raptor=3.0.1"` (already included in ganon install via conda).

!!! Note
    raptor is required to build databases with the Hierarchical Interleaved Bloom Filter (`ganon build --filter-type hibf`)
    To build old style ganon indices `ganon build --filter-type ibf`, raptor is not required

To install raptor from source, follow the instructions below:

#### Dependencies
 
 - CMake >= 3.5
 - GCC 11, 12 or 13 (most recent minor version)

#### Downloading and building raptor + submodules

```bash
git clone --branch raptor-v3.0.1 --recurse-submodules https://github.com/seqan/raptor
```

```bash
cd raptor
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-std=c++23 -Wno-interference-size" ..
make -j 4
```

- binaries will be located in the `bin` directory
- you may have to inform `ganon build` the path to the binaries with `--raptor-path raptor/build/bin`

### Testing

If everything was properly installed, the following command should show the help pages without errors:

```bash
ganon -h
```

#### Running tests

```bash
python3 -m pip install "parameterized>=0.9.0" # Alternative: conda install -c conda-forge "parameterized>=0.9.0"
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/  # optional - downloads large files
cd build/
ctest -VV .
```

## Parameters

```
usage: ganon [-h] [-v]
             {build,build-custom,update,classify,reassign,report,table} ...

- - - - - - - - - -
   _  _  _  _  _   
  (_|(_|| |(_)| |  
   _|   v. 2.1.0
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
usage: ganon build [-h] [-g [...]] [-a [...]] [-l] [-b [...]] [-o] [-c] [-r] [-u] [-m [...]] [-z [...]]
                   [--skip-genome-size] -d DB_PREFIX [-x] [-t] [-p] [-k] [-w] [-s] [-f] [-j] [-y] [-v] [--restart]
                   [--verbose] [--quiet] [--write-info-file]

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

database arguments:
  -l , --level          Highest level to build the database. Options: any available taxonomic rank [species, genus,
                        ...], 'leaves' for taxonomic leaves or 'assembly' for a assembly/strain based analysis (default:
                        species)

download arguments:
  -b [ ...], --source [ ...]
                        Source to download [refseq, genbank] (default: ['refseq'])
  -o , --top            Download limited assemblies for each taxa. 0 for all. (default: 0)
  -c, --complete-genomes
                        Download only sub-set of complete genomes (default: False)
  -r, --representative-genomes
                        Download only sub-set of representative genomes (default: False)
  -u , --genome-updater 
                        Additional genome_updater parameters (https://github.com/pirovc/genome_updater) (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)

important arguments:
  -x , --taxonomy       Set taxonomy to enable taxonomic classification, lca and reports [ncbi, gtdb, skip] (default:
                        ncbi)
  -t , --threads 

advanced arguments:
  -p , --max-fp         Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j , --mode           Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. Only valid for --filter-type ibf. (default:
                        avg)
  -y , --min-length     Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v , --filter-type    Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
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
usage: ganon build-custom [-h] [-i [...]] [-e] [-c] [-n] [-a] [-l] [-m [...]] [-z [...]] [--skip-genome-size] [-r [...]]
                          [-q [...]] -d DB_PREFIX [-x] [-t] [-p] [-k] [-w] [-s] [-f] [-j] [-y] [-v] [--restart]
                          [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). Mutually exclusive --input-file. (default: None)
  -e , --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: fna.gz)
  -c, --input-recursive
                        Look for files recursively in folder(s) provided with --input (default: False)
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (default: None)

custom arguments:
  -n , --input-file     Tab-separated file with all necessary file/sequence information. Fields: file [<tab> target
                        <tab> node <tab> specialization <tab> specialization name]. For details:
                        https://pirovc.github.io/ganon/custom_databases/. Mutually exclusive --input (default: None)
  -a , --input-target   Target to use [file, sequence]. Parse input by file or by sequence. Using 'file' is recommended
                        and will speed-up the building process (default: file)
  -l , --level          Max. level to build the database. By default, --level is the --input-target. Options: any
                        available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy). Further
                        specialization options [assembly, custom]. assembly will retrieve and use the assembly accession
                        and name. custom requires and uses the specialization field in the --input-file. (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)

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
  -p , --max-fp         Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j , --mode           Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. Only valid for --filter-type ibf. (default:
                        avg)
  -y , --min-length     Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v , --filter-type    Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
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
                      [-c [...]] [-e [...]] [-m] [--ranks [...]] [--min-count] [--report-type] [--skip-report] [-o]
                      [--output-one] [--output-all] [--output-unclassified] [--output-single] [-t] [-b] [-f [...]]
                      [-l [...]] [--verbose] [--quiet]

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
                        Min. percentage of a read (set of k-mers) shared with a reference necessary to consider a match.
                        Generally used to remove low similarity matches. Single value or one per database (e.g. 0.7 1
                        0.25). 0 for no cutoff (default: [0.75])
  -e [ ...], --rel-filter [ ...]
                        Additional relative percentage of matches (relative to the best match) to keep. Generally used
                        to keep top matches above cutoff. Single value or one per hierarchy (e.g. 0.1 0). 1 for no
                        filter (default: [0.1])

post-processing/report arguments:
  -m , --multiple-matches 
                        Method to solve reads with multiple matches [em, lca, skip]. em -> expectation maximization
                        algorithm based on unique matches. lca -> lowest common ancestor based on taxonomy. The EM
                        algorithm can be executed later with 'ganon reassign' using the .all file (--output-all).
                        (default: em)
  --ranks [ ...]        Ranks to report taxonomic abundances (.tre). empty will report default ranks [superkingdom,
                        phylum, class, order, family, genus, species, assembly]. (default: [])
  --min-count           Minimum percentage/counts to report an taxa (.tre) [use values between 0-1 for percentage, >1
                        for counts] (default: 5e-05)
  --report-type         Type of report (.tre) [abundance, reads, matches, dist, corr]. More info in 'ganon report'.
                        (default: abundance)
  --skip-report         Disable tree-like report (.tre) at the end of classification. Can be done later with 'ganon
                        report'. (default: False)

output arguments:
  -o , --output-prefix 
                        Output prefix for output (.rep) and tree-like report (.tre). Empty to output to STDOUT (only
                        .rep) (default: None)
  --output-one          Output a file with one match for each read (.one) either an unique match or a result from the EM
                        or a LCA algorithm (--multiple-matches) (default: False)
  --output-all          Output a file with all unique and multiple matches (.all) (default: False)
  --output-unclassified
                        Output a file with unclassified read headers (.unc) (default: False)
  --output-single       When using multiple hierarchical levels, output everything in one file instead of one per
                        hierarchy (default: False)

other arguments:
  -t , --threads        Number of sub-processes/threads to use (default: 1)
  -b, --binning         Optimized parameters for binning (--rel-cutoff 0.25 --rel-filter 0 --min-count 0 --report-type
                        reads). Will report sequence abundances (.tre) instead of tax. abundance. (default: False)
  -f [ ...], --fpr-query [ ...]
                        Max. false positive of a query to accept a match. Applied after --rel-cutoff and --rel-filter.
                        Generally used to remove false positives matches querying a database build with large --max-fp.
                        Single value or one per hierarchy (e.g. 0.1 0). 1 for no filter (default: [1e-05])
  -l [ ...], --hierarchy-labels [ ...]
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
usage: ganon reassign [-h] -i  -o OUTPUT_PREFIX [-e] [-s] [--remove-all] [--skip-one] [--verbose] [--quiet]

options:
  -h, --help            show this help message and exit

required arguments:
  -i , --input-prefix   Input prefix to find files from ganon classify (.all and optionally .rep) (default: None)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix for reassigned file (.one and optionally .rep). In case of multiple files, the
                        base input filename will be appended at the end of the output file 'output_prefix +
                        FILENAME.out' (default: None)

EM arguments:
  -e , --max-iter       Max. number of iterations for the EM algorithm. If 0, will run until convergence (check
                        --threshold) (default: 10)
  -s , --threshold      Convergence threshold limit to stop the EM algorithm. (default: 0)

other arguments:
  --remove-all          Remove input file (.all) after processing. (default: False)
  --skip-one            Do not write output file (.one) after processing. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
```

</details>

<details>
  <summary>ganon report</summary>

```
usage: ganon report [-h] -i [...] [-e INPUT_EXTENSION] -o OUTPUT_PREFIX [-d [...]] [-x] [-m [...]] [-z [...]]
                    [--skip-genome-size] [-f] [-t] [-r [...]] [-s] [-a] [-y] [-p [...]] [-k [...]] [-c] [--verbose]
                    [--quiet] [--min-count] [--max-count] [--names [...]] [--names-with [...]] [--taxids [...]]

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
  --skip-genome-size    Do not attempt to get genome sizes. Valid only without --db-prefix. Activate this option when
                        using sequences not representing full genomes. (default: False)

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
