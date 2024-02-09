# Reports

`ganon report` filters and generates several reports and summaries from the results obtained with `ganon classify`. It is possible to summarize the results in terms of taxonomic and sequence abundances as well as total number of matches.

## Examples

Given the output `.rep` from `ganon classify` and the database used (`--db-prefix`):

### Taxonomic profile with abundance estimation (default)

```bash
ganon report --db-prefix mydb --input results.rep --output-prefix tax_profile --report-type abundance
```

### Sequence profile

```bash
ganon report --db-prefix mydb --input results.rep --output-prefix seq_profile --report-type reads
```

### Matches profile

```bash
ganon report --db-prefix mydb --input results.rep --output-prefix matches --report-type matches
```

### Filtering results

```bash
ganon report --db-prefix mydb --input results.rep --output-prefix filtered --min-count 0.0005 --top-percentile 0.8
```

This will keep only results with a min. abundance of `0.05%` and only the top `80%` most abundant.

## Parameter details

### report type (--report-type)

Several reports are available with `--report-type`: `reads`, `abundance`, `dist`, `corr`, `matches`:

`reads` reports **sequence abundances** which are the basic proportion of reads classified in the sample.

`abundance` will convert sequence abundance into **taxonomic abundances** by re-distributing read counts among leaf nodes and correcting by genome size. The re-distribution applies for reads classified with a LCA assignment and it is proportional to the number of unique matches of leaf nodes available in the ganon database (relative to the LCA node). Genome size is estimated based on [NCBI or GTDB auxiliary files](custom_databases.md#genome-sizes-genome-size-files). Genome size correction is applied by rank based on default ranks only (superkingdom phylum class order family genus species assembly). Read counts in intermediate ranks will be corrected based on the closest parent default rank and re-assigned to its original rank.

`dist` is the same of `reads` with read count re-distribution

`corr` is the same of `reads` with correction by genome size

`matches` will report the total number of matches classified, either unique or shared. *This option will output the total number of matches instead the total number of reads*
