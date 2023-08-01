# Classification

`ganon classify` will match single and/or paired-end sets of reads against one or [more databases](#multiple-and-hierarchical-classification). 
By default, parameters are optimized for **taxonomic profiling**. 

Example:

```bash
ganon classify --db-prefix my_db --paired-reads reads.1.fq.gz reads.2.fq.gz --output-prefix results --threads 32
```

`ganon report` will be automatically executed after `ganon classify` and a [report will be created `.tre`](../outputfiles/#ganon-report).

ganon can perform **taxonomic profiling** and/or **binning** (one tax. assignment for each read) at a taxonomic, strain or sequence level with `ganon classify` + `ganon report`. Some guidelines are listed below, please choose the parameters according to your application:

### Profiling

`ganon classify` is set-up by default to perform taxonomic profiling. It uses:

 - strict thresholds: `--rel-cutoff 0.75` and `--rel-filter 0`

`ganon report` will automatically run after classification with:

 - `--min-count 0.005` (0.5%) to exclude low abundant taxa
 - `--report-type abundance` to generate taxonomic abundances, re-distributing read counts and correcting for genome sizes

!!! Note
    `ganon report` can be used independently from `ganon classify` with the output file `.rep`

### Binning

To achieve better results for taxonomic binning or sequence classification, ganon can be configured with:

 - `--output-all` and `--output-lca` to write `.all` `.lca` files for binning results
 - less strict `--rel-cutoff` and `--rel-filter` values (e.g. `0.25` and `0.1`, respectively)
 - activate the `--reassign` to apply an EM algorithm, re-assigning reads with LCA matches to one most probable target (defined by `--level` in the build procedure). In this case, the `.all` file will be re-generated with one assignment per read.

!!! Note
    `ganon reassign` can be used independently from `ganon classify` with the output file `.rep` and `.all`

!!! tip
    Database parameters can also influence your results. Lower `--max-fp` (e.g. 0.1, 0.001) and higher `--kmer-size` (e.g. `23`, `27`) will improve sensitivity of your results at cost of a larger database and memory usage


## Multiple and Hierarchical classification

`ganon classify` can be performed in multiple databases at the same time. The databases can also be provided in a hierarchical order. 

Multiple database classification can be performed providing several inputs for `--db-prefix`. They are required to be built with the same `--kmer-size` and `--window-size` values. Multiple databases are considered as one (as if built together) and redundancy in content (same reference in two or more databases) is allowed.

To classify reads in a hierarchical order, `--hierarchy-labels` should be provided. When using multiple hierarchical levels, output files will be generated for each level (use `--output-single` to generate a single output from multiple hierarchical levels). Please note that some parameters are set for each database (e.g. `--rel-cutoff`) while others are set for each hierarchical level (e.g. `--rel-filter`)

<details>
  <summary>Examples</summary>
  <br>
Classification against 3 database (as if they were one) using the same cutoff:

```bash
ganon classify --db-prefix db1 db2 db3 \
               --rel-cutoff 0.75 \
               --single-reads reads.fq.gz
```

Classification against 3 database (as if they were one) using different error rates for each:

```bash
ganon classify --db-prefix  db1 db2 db3 \
               --rel-cutoff 0.2 0.3 0.1 \
               --single-reads reads.fq.gz
```

In this example, reads are going to be classified first against db1 and db2. Reads without a valid match will be further classified against db3. `--hierarchy-labels` are strings and are going to be sorted to define the hierarchy order, disregarding input order:

```bash
ganon classify --db-prefix            db1     db2      db3 \
               --hierarchy-labels 1_first 1_first 2_second \
               --single-reads reads.fq.gz
```

In this example, classification will be performed with different `--rel-cutoff` for each database. For each hierarchy levels (`1_first` and `2_second`) a different `--rel-filter` will be used:

```bash
ganon classify --db-prefix            db1     db2      db3 \
               --hierarchy-labels 1_first 1_first 2_second \
               --rel-cutoff             1     0.5     0.25 \
               --rel-filter           0.1              0.5 \
               --single-reads reads.fq.gz
```

</details>
<br>

## Parameter details

### reads (--single-reads, --paired-reads)

ganon accepts single-end and paired-end reads. Both types can be use at the same time. In paired-end mode, reads are always reported with the header of the first pair. Paired-end reads are classified in a standard forward-reverse orientation. The max. read length accepted is 65535 (accounting both reads in paired mode).

### cutoff and filter (--rel-cutoff, --rel-filter)

ganon has two parameters to control a match between reads and references: `--rel-cutoff` and `--rel-filter`.

Every read can be classified against none, one or more references. What will be reported is the remaining matches after `cutoff` and `filter` thresholds are applied, based on the number of shared minimizers (or k-mers) between sequences.

The `cutoff` is the first. It should be set as a minimal value to consider a match between a read and a reference. Next the `filter` is applied to the remaining matches. `filter` thresholds are relative to the best scoring match and control how far from the best match further matches are allowed. `cutoff` can be interpreted as the lower bound to discard spurious matches and `filter` as the fine tuning to control what to keep.

For example, using `--kmer-size 19` (and `--window-size 19` to simplify the example), a certain read (100bp) has the following matches with the 5 references (`ref1..5`), sorted by shared k-mers:

| reference | shared k-mers |
|-----------|---------------|
| ref1      | 82            |
| ref2      | 68            |
| ref3      | 44            |
| ref4      | 25            |
| ref5      | 20            |

this read can have at most 82 shared k-mers (`100-19+1=82`). With `--rel-cutoff 0.25`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.25 |
|-----------|---------------|-------------------|
| ref1      | 82            |                   |
| ref2      | 68            |                   |
| ref3      | 44            |                   |
| ref4      | 25            |                   |
| ~~ref5~~  | ~~20~~        | X                 |

since the `--rel-cutoff` threshold is `82 * 0.25 = 21` (ceiling is applied). Further, with `--rel-filter 0.3`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.25 | --rel-filter 0.3 |
|-----------|---------------|-------------------|------------------|
| ref1      | 82            |                   |                  |
| ref2      | 68            |                   |                  |
| ~~ref3~~  | ~~44~~        |                   | X                |
| ~~ref4~~  | ~~25~~        |                   | X                |
| ~~ref5~~  | ~~20~~        | X                 |                  |


since best match is 82, the filter parameter is removing any match below `0.3 * 82 = 57` (ceiling is applied) shared k-mers. `ref1` and `ref2` are reported as matches.

For databases built with `--window-size`, the relative values are not based on the maximum number of possible shared k-mers but on the actual number of unique minimizers extracted from the read.

A different `cutoff` can be set for every database in a multiple or hierarchical database classification. A different `filter` can be set for every level of a hierarchical database classification.

Note that reads that remain with only one reference match (after `cutoff` and `filter` are applied) are considered a unique match.

### False positive of a query (--fpr-query)

ganon uses Bloom Filters, probabilistic data structures that may return false positive results. The base false positive of a ganon index is controlled by `--max-fp` when building the database. However, this value is the expected false positive for each k-mer. In practice, a sequence (several k-mers) will have a way smaller false positive. ganon calculates the false positive rate of a query as suggested by (Solomon and Kingsford, 2016). The `--fpr-query` will control the max. value accepted to consider a match between a sequence and a reference, avoiding false positives that may be introduce by the properties of the data structure. 

By default, `--fpr-query 1e-5` is used and it is applied after the `--rel-cutoff` and `--rel-filter`. Values between `1e-3` and `1e-10` are recommended. This threshold becomes more important when building smaller databases with higher `--max-fp`, assuring that the false positive is under control. In this case however, you may notice a in sensitivity of your results.

!!! Note
    The false positive of a query was first propose in: Solomon, Brad, and Carl Kingsford. “Fast Search of Thousands of Short-Read Sequencing Experiments.” Nature Biotechnology 34, no. 3 (2016): 1–6. https://doi.org/10.1038/nbt.3442.