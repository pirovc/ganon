
# Important parameters

The most important parameters and trade-offs to be aware of when using ganon:

### ganon build

- `--level`: Highest level to build the database. Can be a taxonomic rank [species, genus, ...], 'leaves' for taxonomic leaves or 'assembly' for a assembly/strain based analysis. The more specific the level, the bigger the database will be.
- `--max-fp`: controls the false positive of the bloom filters. The higher the `--max-fp`, the smaller the databases at a cost of sensitivity in classification.
- `--window-size --kmer-size`: the *window* value should always be the same or larger than the *k-mer* value. The larger the difference between them, the smaller the database will be. However, some sensitivity/precision loss in classification is expected with small *k-mer* and/or large *window*. Larger *k-mer* values (e.g. `31`) will improve classification, specially read binning, at a cost of larger databases.

### ganon classify

- `--rel-cutoff`: defines the min. percentage of k-mers shared to a reference to consider a match. Higher values will improve precision and decrease sensitivity. For taxonomic profiling, a higher value between `0.4` and `0.8` may provide better results. For read binning, lower values between `0.2` and `0.4` are recommended. 
    - **lower** values -> **more read matches**
    - **higher** values -> **less read matches**
- `--rel-filter`: filter matches in relation to the best and worst after the cutoff is applied. `0` means only matches with top score (# of *k-mers*) as the best match will be kept.
    - **lower** values -> **more unique matching reads**
    - **higher** values -> **more multi-matching reads**
- `--multiple-matches`: defines how ganon treats multiple-matching reads. Either by an EM-algorithm based on unique matches or a taxonomy-based LCA algorithm.

### ganon report

!!! tip
    `--report-type` and `--min-count` can be set directly on `ganon classify`

- `--report-type`: reports either taxonomic, sequence or matches abundances. Use `abundance` for taxonomic profiling, `reads` for sequence profiling, or `matches` to report a summary of all matches.
- `--min-count`: cutoff to discard underrepresented taxa. Useful to remove the common long tail of spurious matches and false positives when performing classification. Values between `0.0001` (0.01%) and `0.001` (0.1%) improved sensitivity and precision in our evaluations. The higher the value, the more precise the outcome, with a sensitivity loss. Alternatively `--top-percentile` can be used to keep a relative amount of taxa instead a hard cutoff.

The numeric values above are averages from several experiments with different sample types and database contents. They may not work as expected for your data. If you are not sure which values to use or see something unexpected, please open an [issue](https://github.com/pirovc/ganon/issues){target="_blank"}.