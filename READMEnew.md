# ganon [![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon)[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

ganon classifies short reads against large sets of refence sequences efficiently, with reference download, taxonomic classification (ncbi and gtdb) and many other [features](#Features) [Publication 10.1101/406017](https://dx.doi.org/10.1093/bioinformatics/btaa458)

## Quick install/usage guide

```bash
# Installation from conda
conda install -c bioconda -c conda-forge ganon
# Downloads and builds a database for archaeal complete genome sequences from NCBI RefSeq
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 12
# Classify reads
ganon classify --db-prefix arc_cg_rs --output-prefix classify_results --single-reads my_reads.fq.gz --threads 12
# Generate filtered reports
ganon report --db-prefix arc_cg_rs --input classify_results.rep --output-prefix filtered_report --min-count 0.01
```

## Details

ganon does bla bal bal

### Features

- 1
- 2
- 3

## Installation guide

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

```bash
conda install -c bioconda -c conda-forge ganon
```

* There are possible performance benefits compiling ganon from source in the target machine rather than using the conda version. To do so, please follow the instructions below

<details>
  <summary>Instructions</summary>

### build dependencies

System packages:
- gcc >=8
- cmake >=3.10
- zlib

### run dependencies

System packages:
- python >=3.5
- pandas >=0.22.0
- multitax >=1.1.1

```bash
python3 -m pip install "pandas>=0.22.0"
python3 -m pip install "multitax==1.1.1"
```

### Downloading and building ganon

```bash
git clone --recurse-submodules https://github.com/pirovc/ganon.git
```
  
```bash
cd ganon
python3 setup.py install --record files.txt #optional
mkdir build_cpp
cd build_cpp
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF ..
make
sudo make install #optional
```

- to change install location (e.g. `/myprefix/bin/`), set the installation prefix in the cmake command with `-DCMAKE_INSTALL_PREFIX=/myprefix/ `

- use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.

If everything was properly installed, the following commands should show the help pages without errors:

```bash
ganon -h
```

### Run tests

```bash
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/
cd build_cpp/
ctest -VV .
```

</details>

## Usage guide

save states

### Examples

### Building custom indices

### Updating indices

## Output files

## Multiple and Hierarchical classification

## Choosing/explaining parameters

## Parameters