# Quick Start Guide

## Install

```sh
conda install -c conda-forge -c bioconda ganon
```

## Download and Build a database 

- NCBI RefSeq Bacterial reference genomes

```bash
ganon build --db-prefix bac_rs_rg --source refseq --organism-group bacteria --reference-genomes --threads 24
```

## Classify and generate a tax. profile

- [Download test reads](https://github.com/pirovc/ganon_benchmark/raw/master/files/reads/cami/toy/H01_1M_0.1.fq.gz){target="_blank"}

```bash
ganon classify --db-prefix bac_rs_rg --output-prefix classify_results --single-reads H01_1M_0.1.fq.gz --threads 24
```

- `classify_results.tre` -> taxonomic profile

---
