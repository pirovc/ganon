# ganon [![GitHub release (latest by date)](https://img.shields.io/github/v/release/pirovc/ganon)](https://github.com/pirovc/ganon)

[![Build Status](https://app.travis-ci.com/pirovc/ganon.svg?token=q6Nfx8pLHh8hV3hLz3Pq&branch=master)](https://app.travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1101%2F406017-blue)](https://dx.doi.org/10.1093/bioinformatics/btaa458)

[ganon2 pre-print](https://www.biorxiv.org/content/10.1101/2023.12.07.570547)

ganon2 classifies DNA sequences against large sets of genomic reference sequences efficiently. It features:

- integrated download and build of any subset from RefSeq/Genbank/GTDB with incremental updates
- NCBI and GTDB native support for taxonomic classification, custom taxonomy or no taxonomy at all
- customizable database build for local or non-standard sequence files
- optimized taxonomic binning and classification configurations
- build and classify at various taxonomic levels, strain, assembly, file, sequence or custom specialization
- hierarchical classification using several databases in one or more levels in just one run
- EM and/or LCA algorithms to solve multiple-matching reads
- reporting of multiple and unique matches for every read
- reporting of sequence, taxonomic or multi-match abundances with optional genome size correction
- advanced tree-like reports with several filter options
- generation of contingency tables with several filters for multi-sample studies

Find out more information in the user manual: https://pirovc.github.io/ganon/

## Quick install and usage

```sh
# Install
conda install -c bioconda -c conda-forge ganon
# Download and Build (Archaea - complete genomes - NCBI RefSeq)
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 24
# Classify
ganon classify --db-prefix arc_cg_rs --output-prefix classify_results --paired-reads my_reads.1.fq.gz my_reads.2.fq.gz --threads 24
```

For further examples, database build guides, installation from source and more: https://pirovc.github.io/ganon/