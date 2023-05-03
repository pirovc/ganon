# Classification

## Profiling

## Binning

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