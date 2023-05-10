# ganon [![GitHub release (latest by date)](https://img.shields.io/github/v/release/pirovc/ganon)](https://github.com/pirovc/ganon) [![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1101%2F406017-blue)](https://dx.doi.org/10.1093/bioinformatics/btaa458)

ganon classifies DNA sequences against large sets of genomic reference sequences efficiently. It features:

- automatic download, build and update procedures for commonly used databases (RefSeq and GenBank)
- classification with binning and taxonomic profiling
- multiple taxonomy integration (NCBI and GTDB) with lowest common ancestor (LCA)
- read reassignment EM algorithm for multi-matching reads
- hierarchical use of multiple databases
- taxonomic and sequence abundance reports with genome size correction
- advanced reporting and filtration of results
- contingency table generation

Find out more information in the user manual: https://pirovc.github.io/ganon/

## Quick install with conda

```sh
conda install -c bioconda -c conda-forge ganon
```

## Basic usage

### Download and Build (Archaea - complete genomes - NCBI RefSeq)

```bash
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 24
```

### Classify
```bash
ganon classify --db-prefix arc_cg_rs --output-prefix classify_results --paired-reads my_reads.1.fq.gz my_reads.2.fq.gz --threads 24
```

For further examples, database guide, installation from source and more: https://pirovc.github.io/ganon/