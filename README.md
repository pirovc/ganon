# ganon2 

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/pirovc/ganon)](https://github.com/pirovc/ganon) [![Build Status](https://app.travis-ci.com/pirovc/ganon.svg?token=q6Nfx8pLHh8hV3hLz3Pq&branch=main)](https://app.travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/main/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1093/nargab/lqaf094-blue)](https://dx.doi.org/10.1093/nargab/lqaf094)

ganon2 index large sets of genomic reference sequences efficiently and quickly classify sequence fragments to their closest matching reference, providing taxonomic or sequence abundance estimations.

Find out more information in the user manual: https://pirovc.github.io/ganon/

## Quick install and usage

```sh
# Install
conda install -c conda-forge -c bioconda ganon
# Download and Build (Bacteria - complete genomes - NCBI RefSeq)
ganon build --db-prefix bac_cg_rs --source refseq --organism-group bacteria --complete-genomes --threads 24 --download-threads 12
# Classify
ganon classify --db-prefix bac_cg_rs --output-prefix classify_results --paired-reads my_reads.1.fq.gz my_reads.2.fq.gz --threads 24
```

For further examples, database build guides, installation from source and more: https://pirovc.github.io/ganon/