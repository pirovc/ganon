# ganon [![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1101%2F406017-blue)](https://dx.doi.org/10.1093/bioinformatics/btaa458)

ganon classifies DNA sequences against large sets of genomic reference sequences efficiently. It features:

- automatic download, build and update for commonly used repos (refseq/genbank)
- binning and taxonomic profiling
- multiple taxonomy integration (ncbi/gtdb)
- LCA algorithm + read reassignment EM algorithm
- hierarchical use of databases
- taxonomic and sequence abundance reports with genome size correction
- contingency tables and many more

Documentation: https://pirovc.github.io/ganon/

## Quick install

```sh
conda install -c bioconda -c conda-forge ganon
```

## Basic usage

### Download and build
```bash
# Archaeal complete genome sequences from NCBI RefSeq
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 12
```

### Classify
```bash
ganon classify --db-prefix arc_cg_rs --output-prefix classify_results --single-reads my_reads.fq.gz --threads 12
```

### Re-generate reports and create tables from multiple reports
```bash
ganon report --db-prefix arc_cg_rs --input classify_results.rep --output-prefix filtered_report --min-count 0.01
ganon table --input classify_results.tre filtered_report.tre --output-file output_table.tsv --top-sample 10
```

### Update the database at a later time point
```bash
ganon update --db-prefix arc_cg_rs --threads 12
```
