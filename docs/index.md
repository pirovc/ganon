# ganon2

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/pirovc/ganon)](https://github.com/pirovc/ganon){target="_blank"} [![Build Status](https://app.travis-ci.com/pirovc/ganon.svg?token=q6Nfx8pLHh8hV3hLz3Pq&branch=main)](https://app.travis-ci.com/pirovc/ganon){target="_blank"} [![codecov](https://codecov.io/gh/pirovc/ganon/branch/main/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon){target="_blank"} [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon){target="_blank"} [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon){target="_blank"} [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/ganon/README.html){target="_blank"} [![Publication](https://img.shields.io/badge/DOI-10.1093/nargab/lqaf094-blue)](https://dx.doi.org/10.1093/nargab/lqaf094){target="_blank"}

Code: [GitHub repository](https://github.com/pirovc/ganon){target="_blank"}

ganon index large sets of genomic reference sequences efficiently and quickly classify sequence fragments to their closest matching reference. It uses [Hierarchical Interleaved Bloom Filters](https://doi.org/10.1186/s13059-023-02971-4){target="_blank"} as indices based on k-mers/minimizers. After classification, taxonomic or sequence abundances are estimated and reported.

## Features

- [Default database build](default_databases.md/)
    - integrated download of any subset from [RefSeq/Genbank/GTDB](default_databases.md#refseq-and-genbank)
    - incremental [updates](default_databases.md#update-ganon-update)
    - [multi-taxonomy](custom_databases.md/#taxonomy) support: [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy){target="_blank"}, [GTDB](default_databases.md#gtdb){target="_blank"}
- [Custom database build](custom_databases.md)
    - using local or [non-standard](custom_databases.md/#non-standard-filesheaders-with-input-file) sequence files
    - supports all taxonomic [levels](custom_databases.md/#level) plus strain, assembly, file, sequence or custom specialization
    - [conversion](custom_databases.md/#conversion) between NCBI<>GTDB and GTDB versions
- [Classification](classification.md/) 
    - [taxonomic binning](classification.md#binning) and [profiling](classification.md#profiling) configurations
    - [batch classification](classification.md/#batch-classification) processing several samples in one run
    - [hierarchical classification](classification.md#multiple-and-hierarchical-classification) to classify against several databases in one or more levels
    - [EM and/or LCA](classification.md#reads-with-multiple-matches) algorithms to solve multiple-matching reads
- [Reports](reports.md/) 
    - reports of multiple and unique matches for every read
    - [reporting](reports.md#report-types-report-type) of sequence, taxonomic or match abundances with optional genome size correction
    - tree-like [reports](reports.md) with several filter options
- [Contingency tables](table.md) with several filters for multi-sample studies

ganon achieved very good results in [our own evaluations](https://dx.doi.org/10.1093/bioinformatics/btaa458){target="_blank"} but also in independent evaluations: [LEMMI](https://lemmi-v1.ezlab.org/){target="_blank"}, [LEMMI v2](https://lemmi.ezlab.org/){target="_blank"} and [CAMI2](https://dx.doi.org/10.1038/s41592-022-01431-4){target="_blank"} (taxonomic [profiling](https://cami-challenge.org/taxonomic_profiling/){target="_blank"} and [binning](https://cami-challenge.org/taxonomic_binning/){target="_blank"}).

## Installation with conda

The easiest way to install ganon is via conda, using the bioconda and conda-forge channels:

```bash
conda install -c conda-forge -c bioconda ganon
```

## Installation from source

### Dependencies

#### Python

- python >=3.10
- pandas >=1.2.0
- [multitax](https://github.com/pirovc/multitax){target="_blank"} >=1.5.1
- [genome_updater](https://github.com/pirovc/genome_updater){target="_blank"} >=0.7.2

#### C++

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
cd ganon
pip install .

# Compile and install C++ side
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF -DLONGREADS=OFF ..
make -j 4
sudo make install  # optional, otherwise indicate path when running ganon with --ganon-path
```

- to change install location (e.g. `/myprefix/bin/`), set the installation prefix in the cmake command with `-DCMAKE_INSTALL_PREFIX=/myprefix/ `
- use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.
- to classify extremely large reads or contigs that would need more than 65000 k-mers, use `-DLONGREADS=ON`

### Installing raptor

The easiest way to install [raptor](https://github.com/seqan/raptor){target="_blank"} is via conda with `conda install -c bioconda -c conda-forge "raptor=3.0.1"` (already included in ganon install via conda).

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
pip install .[dev]
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/  # optional - downloads large files
cd build/
ctest -VV .
```