# Reports

`ganon report` generates taxonomic reports and summaries based on the results obtained using `ganon classify`. Optionally filters and formats the report. Results results can be summarised in terms of taxonomic and sequence abundances, as well as the total number of matches.

A file with `.tre` extension is generated containing the taxonomic report, more infos about the format [here](outputfiles.md/#ganon-report). Additionaly, a summary is reported (STDERR) ny taxonomic rank, for example:

```txt
         unique     shared     children   total      
root     0%         0%         71.75%     71.75%     
domain   0%         0%         71.75%     71.75%     
phylum   0%         0%         71.75%     71.75%     
class    0%         0%         71.55%     71.55%     
order    0%         0%         71.7%      71.7%      
family   0%         0%         71.6%      71.6%      
genus    0%         0%         71.75%     71.75%     
species  26.65%     16.25%     28.85%     71.75%     
assembly 17.2%      11.65%     0%         28.85%  
```

Details on unique, shared and children values can be found [here](outputfiles.md/#ganon-report).

## Report types (--report-type)

Several reports are available with `--report-type`: `reads`, `abundance`, `dist`, `corr`, `matches`:

- `reads` reports **sequence abundances** which are the basic proportion of reads classified in the sample.
- `abundance` will convert sequence abundance into **taxonomic abundances** by re-distributing read counts among leaf nodes and correcting by genome size. The re-distribution applies for reads classified with a LCA assignment and it is proportional to the number of unique matches of leaf nodes available in the ganon database (relative to the LCA node). If EM was previously used to re-distribute reads, will only correct for sequence abundance (same as `corr`). Genome size is estimated based on [NCBI or GTDB auxiliary files](custom_databases.md#genome-sizes-genome-size-files). Genome size correction is applied by rank based on default ranks only (superkingdom phylum class order family genus species assembly). Read counts in intermediate ranks will be corrected based on the closest parent default rank and re-assigned to its original rank.
- `dist` is the same of `reads` with read count re-distribution
- `corr` is the same of `reads` with correction by genome size
- `matches` will report the total number of matches classified, either unique or shared. *This option will output the total number of matches instead the total number of reads*

## Filter and format

Results can be filtered in several ways:

- `--ranks` filters specific taxonomic ranks.
- `--top-percentile` keeps only the top percentile chosen for each rank (e.g. `0.5` will keep 50% of taxa for each rank).
- `--min-count/--max-count` min/max percentage or number of counts to keep an taxa. Values between 0-1 for percentage, >1 for number of counts (total).
- `--names/--names-with/--taxids` keeps only entries matching specific names or taxids.
- `--skip-hierarchy/--keep-hierarchy` skip or keep hiearchy levels on the report

Results can be formatted with:

- `--output-format` file output format: `text, tsv, csv, bioboxes`.
- `--sort` order of the output: `rank, lineage, count, unique`.
- `--normalize` will ignore the number of unclassified reads, normalizing the output to 100%.
- `--split-hierarchy` splits output reports by hierarchy levels.

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

