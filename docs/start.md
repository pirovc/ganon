# Quick Start Guide

## Install

```sh
conda install -c bioconda -c conda-forge ganon
```

## Download and Build a database 

- Bacteria - NCBI RefSeq - representative genomes

```bash
ganon build --db-prefix bac_rs_rg --source refseq --organism-group bacteria --representative-genomes --threads 24
```

- If you want to test ganon functionalities with a smaller database, use `archaea` instead of `bacteria` in the example above.

## Classify and generate a tax. profile

- [Download test reads](https://github.com/pirovc/ganon_benchmark/raw/master/files/reads/cami/toy/H01_1M_0.1.fq.gz)

```bash
ganon classify --db-prefix bac_rs_rg --output-prefix classify_results --single-reads H01_1M_0.1.fq.gz --threads 24
```

- `classify_results.tre` -> taxonomic profile

---


## Important parameters

The most important parameters and trade-offs to be aware of when using ganon:

### ganon build

- `--max-fp --filter-size`: controls the false positive of the bloom filters and the size of the filter (which is the same as the amount of memory needed). The higher the `--max-fp`, the smaller the databases at a cost of sensitivity in classification. `--filter-size` can be used instead of `--max-fp` to define a specific size for your database. In this case, the false positive will be reported at the end of the build.
- `--window-size --kmer-size`: the *window* value should always be the same or larger than the *k-mer* value. The larger the difference between them, the smaller the database will be. However, some sensitivity/precision loss in classification is expected with small *k-mer* and/or large *window*. Larger *k-mer* values (e.g. `31`) will improve classification, specially read binning, at a cost of larger databases.
- `--hibf`: build smaller databases that can be queried faster. Building will take longer.

### ganon classify

- `--rel-cutoff`: defines the min. percentage of k-mers shared to a reference to consider a match. Higher values will improve precision and decrease sensitivity. For taxonomic profiling, a higher value between `0.4` and `0.8` may provide better results. For read binning, lower values between `0.2` and `0.4` are recommended. 
    - **lower** values -> **more read matches**
    - **higher** values -> **less read matches**
- `--rel-filter`: filter matches in relation to the best and worst after the cutoff is applied. `0` means only matches with top score (# of *k-mers*) as the best match will be kept.
    - **lower** values -> **more unique matching reads**
    - **higher** values -> **more multi-matching reads**
- `--reassign`: runs an EM-algorithm to reassign reads that received multiple matches. It provides a unique match for each read at the level the database was built (e.g. assembly or species). Mostly useful for read binning, with little overall impact on taxonomic profiling. Can be used independently with `ganon reassign`.

### ganon report

- `--report-type`: reports either taxonomic, sequence or matches abundances. Use `corr` or `abundance` for taxonomic profiling, `reads` or `dist` for sequence profiling and `matches` to report a summary of all matches.
- `--min-count`: cutoff to discard underrepresented taxa. Useful to remove the common long tail of spurious matches and false positives when performing classification. Values between `0.0001` (0.01%) and `0.001` (0.1%) improved sensitivity and precision in our evaluations. The higher the value, the more precise the outcome, with a sensitivity loss. Alternatively `--top-percentile` can be used to keep a relative amount of taxa instead a hard cutoff.

The numeric values above are averages from several experiments with different sample types and database contents. They may not work as expected for your data. If you are not sure which values to use or see something unexpected, please open an [issue](https://github.com/pirovc/ganon/issues).