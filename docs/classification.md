# Classification

`ganon classify` will match single and/or paired-end sets of reads against one or [more databases](#multiple-and-hierarchical-classification). 
By default, parameters are optimized for **taxonomic profiling**, meaning that less reads will be classified but with a higher sensitivity. For example:

```bash
ganon classify --db-prefix my_db --paired-reads reads.1.fq.gz reads.2.fq.gz --output-prefix results --threads 32
```

Output files:

 - `results.rep`: plain report of the run, used to further generate tree-like reports
 - `results.tre`: tree-like report with cumulative abundances by taxonomic ranks (can be re-generated with `ganon report`)

By default, `ganon classify` only write report files. To get files with the classification of each read, use `--output-one` and/or `--output-all`. More information about output files [here](outputfiles.md#ganon-classify).

!!! Note
    ganon performs **taxonomic profiling** and/or **binning** (one tax. assignment for each read) at a taxonomic, strain or sequence level. Some guidelines are listed below, please choose the parameters according to your application.

### Profiling

`ganon classify` is set-up by default to perform taxonomic profiling. It uses:

 - strict thresholds: `--rel-cutoff 0.75` and `--rel-filter 0.1`
 - `--min-count 0.00005` (0.005%) to exclude very low abundant taxa
 - `--report-type abundance` to generate taxonomic abundances, correcting for genome sizes  (more infos [here](reports.md#report-type-report-type))

### Binning

To achieve better results for taxonomic binning or sequence classification, `ganon classify` can be configured with `--binning`, that is the same as:

 - less strict thresholds: `--rel-cutoff 0.25 --rel-filter 0`
 - `--min-count 0` reports all taxa with at least one read assigned to it
 - `--report-type reads` will report sequence abundances instead of taxonomic abundances (more infos [here](reports.md#report-type-report-type))

!!! Tip
    Database parameters in `ganon build` can also influence your results. Lower `--max-fp` (e.g. 0.1, 0.001) and higher `--kmer-size` (e.g. `23`, `27`) will improve sensitivity of your results at cost of a larger database and memory usage.

## Reads with multiple matches

There are two ways to solve reads with multiple-matches in `ganon classify`:

 - `--multiple-matches em` (default): uses an Expectation-Maximization algorithm, re-assigning reads with multiple matches to one most probable target (defined by `--level` in the build procedure).
 - `--multiple-matches lca`: uses the Lowest Common Ancestor algorithm, re-assigning reads with multiple matches to higher common ancestors in the taxonomic tree.
 - `--multiple-matches skip`: will not resolve multi-matching reads

!!! Tip
    - The Expectation-Maximization can be performed independently with `ganon reassign` using the output files `.rep` and `.all`.
    - Reports can be generated independently with `ganon report` using the output file `.rep`

!!! Note
    `--multiple-matches lca` paired with `--report-type abundance` or `dist` will distribute read **counts** with multiple matches to one most probable target (defined by `--level` in the build procedure), instead of a higher taxonomic rank. In this case the distribution is simply based on the number of taxa with unique matches and it is not as precise as the EM algorithm, but it will run faster since the per-read basis re-assignment can be skipped.

## Classifying more reads

By default ganon will classify less reads in favour of sensitivity. To classify more reads, use less strict `--rel-cutoff` and `--rel-filter` values (e.g. `0.25` and `0`, respectively). More details [here](#cutoff-and-filter-rel-cutoff-rel-filter).

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

ganon has two main parameters to control the number and strictess of matches between reads and references: `--rel-cutoff` and `--rel-filter`.

Every read can be classified against none, one or more references. ganon will report only matches after `cutoff` and `filter` thresholds are applied, based on the number of shared k-mers between sequences (use `--rel-cutoff 0` and `--rel-filter 1` to deactivate them).

The `cutoff` is the first to be applied. It sets the min. percentage of k-mers of a read to be shared with a reference to consider a match. Next the `filter` is applied to the remaining matches. `filter` thresholds are relative to the best and worst scoring match after `cutoff` and control the percentage of additional matches (if any) should be reported, sorted from the best to worst. `filter` won't change the total number of matched reads but will change the amount of unique or multi-matched reads.

`cutoff` can be interpreted as the lower bound to discard spurious matches and `filter` as the fine tuning to control what to keep.

In summary:

  - `--rel-cutoff` controls the strictness of the matching algorithm.
    - **lower** values -> **more read matches**
    - **higher** values -> **less read matches**
  - `--rel-filter` controls how many matches each read will have, from best to worst
    - **lower** values -> **more unique matching reads**
    - **higher** values -> **more multi-matching reads**

For example, using a hypothetical number of k-mer matches, a certain read with 82 k-mers has the following matches with the 5 references (`ref1..5`), sorted number of shared k-mers:

| reference | shared k-mers |
|-----------|---------------|
| ref1      | 82            |
| ref2      | 68            |
| ref3      | 44            |
| ref4      | 25            |
| ref5      | 20            |

With `--rel-cutoff 0.25`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.25 |
|-----------|---------------|-------------------|
| ref1      | 82            |                   |
| ref2      | 68            |                   |
| ref3      | 44            |                   |
| ref4      | 25            |                   |
| ~~ref5~~  | ~~20~~        | X                 |

since the `--rel-cutoff` threshold is `82 * 0.25 = 21` (ceiling is applied).

Next, with `--rel-filter 0.5`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.25 | --rel-filter 0.5 |
|-----------|---------------|-------------------|------------------|
| ref1      | 82            |                   |                  |
| ref2      | 68            |                   |                  |
| ~~ref3~~  | ~~44~~        |                   | X                |
| ~~ref4~~  | ~~25~~        |                   | X                |
| ~~ref5~~  | ~~20~~        | X                 |                  |


since `82` is the best match and `25` is the worst remaining match, the filter will keep the top the remaining matches, based on the shared k-mers threshold `82 - ((82-25)*0.5) = 54` (ceiling is applied). `ref1` and `ref2` are reported as matches


!!! Tip
    The actual number of unique k-mers in a read are used as an upper bound to calculate the thresholds. The same is applied when using `--window-size` and minimizers.

!!! Note
    A different `--rel-cutoff` can be set for every database in a multiple or hierarchical database classification. A different `--rel-filter` can be set for every level of a hierarchical database classification.

!!! Note
    Reads that remain with only one reference match (after `cutoff` and `filter` are applied) are considered a unique match.



### False positive of a query (--fpr-query)

ganon uses Bloom Filters, probabilistic data structures that may return false positive results. The base false positive of a ganon index is controlled by `--max-fp` when building the database. However, this value is the expected false positive for each k-mer. In practice, a sequence (several k-mers) will have a way smaller false positive. ganon calculates the false positive rate of a query as suggested by (Solomon and Kingsford, 2016). The `--fpr-query` will control the max. value accepted to consider a match between a sequence and a reference, avoiding false positives that may be introduce by the properties of the data structure. 

By default, `--fpr-query 1e-5` is used and it is applied after the `--rel-cutoff` and `--rel-filter`. Values between `1e-3` and `1e-10` are recommended. This threshold becomes more important when building smaller databases with higher `--max-fp`, assuring that the false positive is under control. In this case however, sensitivity of results may decrease.

!!! Note
    The false positive of a query was first propose in: Solomon, Brad, and Carl Kingsford. “Fast Search of Thousands of Short-Read Sequencing Experiments.” Nature Biotechnology 34, no. 3 (2016): 1–6. https://doi.org/10.1038/nbt.3442.