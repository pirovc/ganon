# ganon [![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon)

a k-mer based DNA classification tool using Interleaved Bloom Filters for short reads

> **ganon: precise metagenomics classification against large and up-to-date sets of reference sequences**
> Vitor C. Piro, Temesgen H. Dadi, Enrico Seiler, Knut Reinert, and Bernhard Y. Renard
> Bioinformatics 2020 doi: [10.1101/406017](https://dx.doi.org/10.1093/bioinformatics/btaa458)

## Details

ganon is a tool designed to index large sets of genomic reference sequences and to classify short reads against them efficiently. Reads are classified and filtered based on shared k-mers or minimizers with a defined absolute error rate using the k-mer counting lemma (q-gram lemma) or a relative cutoff. The tool was mainly developed, but not limited, to the metagenomic classification problem: assign short fragments to their closest reference among thousands of references.

ganon can further:
 - [update indices](#updating-the-index) incrementally in a fraction of time used to built them
 - use [minimizers](#minimizers-or-offset) to build small indices and fast classification 
 - perform [hierarchical classification](#hierarchical-classification)
 - classify at strain or taxonomic level but also at a custom [specialization](#--specialization)
 - integrated NCBI taxonomy to build and classify
 - report the lowest common ancestor (LCA) and/or multiple matches
 - report unique matches between reads and references
 - generate reports and tables for multi-sample studies

ganon relies mainly on 3 external tools/libraries
 - [SeqAn3](https://github.com/seqan/seqan3) and the Interleaved Bloom Filter implementation
 - [TaxSBP](https://github.com/pirovc/taxsbp) to cluster sequences into taxonomic related groups
 - [genome_updater](https://github.com/pirovc/genome_updater) (optional) to download reference sequences 

ganon achieved very good results in [our evaluations](https://dx.doi.org/10.1093/bioinformatics/btaa458) but also in independent evaluations: [LEMMI](https://lemmi.ezlab.org/) and [CAMI2](https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

```bash
conda install -c bioconda -c conda-forge ganon
```

* There are possible performance benefits compiling ganon from source rather than using the conda version. To do so, please follow the instructions: [Install without conda](#install-without-conda)

## Running ganon with sample data

### build

Create index/database from reference genomic sequences

```bash
ganon build --db-prefix sample_bacteria \
            --input-files tests/ganon/data/build/bacteria_*.fasta.gz \
            --taxdump-file tests/ganon/data/mini_nodes.dmp tests/ganon/data/mini_names.dmp \
            --seq-info-file tests/ganon/data/build/bacteria_seqinfo.txt
```

- `ganon build` (with a space in between) should be used and it is different from the `ganon-build` command
- `--taxdump-file` and `--seq-info-file` files are optional and will be automatically downloaded/generated if not provided
- To build a smaller filter that will enable very fast classification (at a small cost of sens./prec.), use minimizers:

<details>
  <summary>Example</summary>

```bash
ganon build --db-prefix sample_bacteria \
            --input-files tests/ganon/data/build/bacteria_*.fasta.gz \
            --taxdump-file tests/ganon/data/mini_nodes.dmp tests/ganon/data/mini_names.dmp \
            --seq-info-file tests/ganon/data/build/bacteria_seqinfo.txt \
            --window-size 32
```

  </details>

### classify

Classify reads against built index/database

```bash
ganon classify --db-prefix sample_bacteria --single-reads tests/ganon/data/bac.sim.1.fq -o sample_results
```

### report

Generate readable reports, with filter options

```bash
ganon report --db-prefix sample_bacteria --rep-file sample_results.rep \
             --min-percentage 0.5 --ranks all \
             --output-prefix new_report -f tsv
```

### table

Generate a rank table (e.g.: samples X species) from multiple reports, with filter options

```bash
ganon table -i sample_results.tre new_report.tre -o sample_table.tsv
```

### update

Update index/database adding and/or removing reference genomic sequences

```bash
ganon update --db-prefix sample_bacteria --output-db-prefix sample_bacteria_virus --input-files tests/ganon/data/update/virus_*.fasta.gz
```

## Building custom indices

To build custom indices, ganon requires one (or multiple) fasta file(s) with the standard NCBI accession.version header (e.g. `>NC_009515.1`). For every sequence, taxonomic information will be automatically retrieved. Check the parameter `--rank`, `--specialization` and `--seq-info`  for more information.

If you are working with sequences from a different source than NCBI, please provide a `--seq-info-file`.

### Downloading sequences and building a new index

We suggest using [genome_updater](https://github.com/pirovc/genome_updater) (`conda install -c bioconda genome_updater`) to download sequences from RefSeq/Genbank. genome_updater can download and keep a subset of those repositories updated. For example:

Downloading archaeal and bacterial complete genomes from RefSeq: 

```bash
genome_updater.sh -g "archaea,bacteria" \
                  -d "refseq" \
                  -l "Complete Genome" \
                  -f "genomic.fna.gz,assembly_report.txt" \
                  -o "RefSeqCG_arc_bac" -b "v1" \
                  -a -m -u -r -p -t 24
```
Where `-a` will additionally download the taxdump.tar.gz, `-m` will force the MD5 check for every file downloaded, `-u -r -p` will generate reports which can be used later, `-t` set the number of parallel downloads, `-b` will name the version and `-o` set the output working directory. It also possible to download references for specific taxonomic groups with genome_updater.

Building the index based on the files of the example above:

```bash
ganon build --db-prefix ganon_db \
            --input-files RefSeqCG_arc_bac/v1/files/*genomic.fna.gz
```

- If you are getting the bash error `Argument list too long` use `--input-directory "RefSeqCG_arc_bac/v1/files/" --input-extension "genomic.fna.gz"` instead of `--input-files`.
- To speed up the building process, consider [using extra files to build and update](#using-extra-files-to-build-and-update)

### Updating the index

To update the folder with the most recent releases (after some days), just re-run the same download command:

```bash
genome_updater.sh -g "archaea,bacteria" \
                  -d "refseq" \
                  -l "Complete Genome" \
                  -f "genomic.fna.gz,assembly_report.txt" \
                  -o "RefSeqCG_arc_bac" -b "v2" \
                  -a -m -u -r -p -t 24
```

This is now going to download the new files (and skip the outdated ones) into the new label `v2`. New log and report files with the current version will be generated. It also checks if the last downloaded version is complete and downloads any missing file.

Updating the index based on the files of the example above:

```bash
ganon update --db-prefix ganon_db --output-db-prefix ganon_db_updated \
             --input-files $(find RefSeqCG_arc_bac/v2/files/ -name *genomic.fna.gz -type f)
```

If `--output-db-prefix` is not set, the database files `ganon_db.*` will be overwritten with the updated version. The `find` command will only look for files `-type f` inside the `v2/files/`, ignoring symbolic links from sequences which were not changing in this version.

By default, `ganon update` will only add new sequences provided to the index. To perform a full update, also removing sequences for the index, use the option `--update-complete` and provide the full set of updated references instead of only the new sequences, as below:

```bash
ganon update --db-prefix ganon_db --output-db-prefix ganon_db_updated \
             --update-complete \
             --input-files RefSeqCG_arc_bac/v2/files/*genomic.fna.gz
```

### Using extra files to build and update

<details>
  <summary>Example</summary>

Optionally, some extra files generated by genome_updater can be further used to speed-up the building process:

```bash
# Extract taxonomic information from genome_updater reports
awk 'BEGIN {FS="\t";OFS="\t"}{if($4!="na"){ print $4,$5,$6,$2 }}' RefSeqCG_arc_bac/v1/updated_sequence_accession.txt > RefSeqCG_arc_bac/v1/seqinfo.txt

# Use generated files from genome_updater on ganon build
ganon build --db-prefix ganon_db \
            --input-files RefSeqCG_arc_bac/v1/files/*genomic.fna.gz \
            --seq-info-file RefSeqCG_arc_bac/v1/seqinfo.txt \
            --taxdump-file RefSeqCG_arc_bac/v1/{TIMESTAMP}_taxdump.tar.gz
```

The same goes for the update process:

```bash
# Extract taxonomic information from genome_updater reports
awk 'BEGIN {FS="\t";OFS="\t"}{if($1=="A" && $4!="na"){ print $4,$5,$6,$2 }}' RefSeqCG_arc_bac/v2/updated_sequence_accession.txt > RefSeqCG_arc_bac/v2/seqinfo.txt

# Use generated files on ganon update
ganon update --db-prefix ganon_db \
             --output-db-prefix ganon_db_updated \
             --input-files $(find RefSeqCG_arc_bac/v2/files/ -name *genomic.fna.gz -type f) \
             --seq-info-file RefSeqCG_arc_bac/v2/seqinfo.txt \
             --taxdump-file RefSeqCG_arc_bac/v2/{TIMESTAMP}_taxdump.tar.gz
```

Obs:
-  If `-d genbank` was used in genome_updater, change the occurrences of `$4` to `$3` in the `awk` commands above.
- `{TIMESTAMP}` should be the timestamp (YYYY-MM-DD_HH-MM-SS) automatically generated when the files were downloaded.

</details>

## Output files

### build/update

Every run on `ganon build` or `ganon update` will generate the following database files:

 - {prefix}**.ibf**: main interleaved bloom filter file
 - {prefix}**.map**: tab-separated mapping between targets and bin identifiers. Targets should be present in the .tax file as a node *(fields: target, bin id)*
 - {prefix}**.tax**: taxonomic tree *(fields: target/node, parent, rank, name)*
 - {prefix}**.gnn**: gzipped pickled file (python) with information about clustering and parameters used

Obs:
-  Database files generated with version 1.0.0 or higher are not compatible with older versions.

### classify
 
 - {prefix}**.rep**: plain report of the run with only targets that received a match *(fields: 1) hierarchy_label, 2) target, 3) total matches, 4) unique reads, 5) lca reads, 6) rank, 7) name)*. At the end prints 2 extra lines with `#total_classified` and `#total_unclassified`
 - {prefix}**.lca**: output with one match for each classified read after LCA. Only generated with `--output-lca` active. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierarchy}.lca *(fields: read identifier, target, (max) k-mer count)*
 - {prefix}**.all**: output with all matches for each read. Only generated with `--output-all` active **Warning: file can be very large**. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierarchy}.all *(fields: 1) read identifier, 2) target, 3) k-mer count)*
  - {prefix}**.tre**: report file (see below)

### report

 - {prefix}**.tre**: tab-separated tree-like report with cumulative counts and taxonomic lineage. By default, this is a read-based report where each read classified is counted once. It is possible to generate this for all read matches (`ganon report --report-type matches`). In this case, single and shared matches are reported to their target. Each line in this report is a taxonomic entry, with the following fields: 

  1) taxonomic rank *(e.g. phylum, species, ...)*
  2) target *(e.g. taxid/specialization)*
  3) target lineage *(e.g 1|2|1224|...)*
  4) target name *(e.g. Paenibacillus polymyxa)*
  5) \# unique assignments *(number of reads that matched exclusively to this target)*
  6) \# shared assignments *(number of reads with non-unique matches directly assigned to this target. Represents the lca matches (`--report-type reads`) or shared matches (`--report-type matches`))*
  7) \# children assignments *(number of reads assigned to all children nodes of this target)*
  8) \# cumulative assignments *(the sum of the unique, shared and children reads/matches assigned up-to this target)*
  9) \% cumulative assignments

- Using `--report-type reads` the first line of the file will show the number of unclassified reads

- The sum of cumulative assignments for the unclassified and root lines should be 100%. The final cumulative sum of reads/matches may be under 100% if any filter is successfully applied and/or hierarchical selection is selected (keep/skip/split).

- When `--report-type reads` only taxa that received direct read matches, either unique or through lca, are considered. Some reads may have only shared matches and will not be reported directly (but will be accounted on some parent level). To look at those matches you can create a report with `--report-type matches` or look at the file {prefix}**.rep**.

### table

 - {output_file}: a tab-separated file with counts/percentages of taxa for multiple samples
 
### Examples of output files

The main output file is the `{prefix}.tre` which will summarize the results:

```
unclassified                                                 unclassified             0   0  0   2   2.02020
root          1       1                                      root                     0   0  97  97  97.97980
superkingdom  2       1|2                                    Bacteria                 0   0  97  97  97.97980
phylum        1239    1|2|1239                               Firmicutes               0   0  57  57  57.57576
phylum        1224    1|2|1224                               Proteobacteria           0   0  40  40  40.40404
class         91061   1|2|1239|91061                         Bacilli                  0   0  57  57  57.57576
class         28211   1|2|1224|28211                         Alphaproteobacteria      0   0  28  28  28.28283
class         1236    1|2|1224|1236                          Gammaproteobacteria      0   0  12  12  12.12121
order         1385    1|2|1239|91061|1385                    Bacillales               0   0  57  57  57.57576
order         204458  1|2|1224|28211|204458                  Caulobacterales          0   0  28  28  28.28283
order         72274   1|2|1224|1236|72274                    Pseudomonadales          0   0  12  12  12.12121
family        186822  1|2|1239|91061|1385|186822             Paenibacillaceae         0   0  57  57  57.57576
family        76892   1|2|1224|28211|204458|76892            Caulobacteraceae         0   0  28  28  28.28283
family        468     1|2|1224|1236|72274|468                Moraxellaceae            0   0  12  12  12.12121
genus         44249   1|2|1239|91061|1385|186822|44249       Paenibacillus            0   0  57  57  57.57576
genus         75      1|2|1224|28211|204458|76892|75         Caulobacter              0   0  28  28  28.28283
genus         469     1|2|1224|1236|72274|468|469            Acinetobacter            0   0  12  12  12.12121
species       1406    1|2|1239|91061|1385|186822|44249|1406  Paenibacillus polymyxa   57  0  0   57  57.57576
species       366602  1|2|1224|28211|204458|76892|75|366602  Caulobacter sp. K31      28  0  0   28  28.28283
species       470     1|2|1224|1236|72274|468|469|470        Acinetobacter baumannii  12  0  0   12  12.12121
```

running `ganon classify` or `ganon report` with `--ranks all`, the output will show all ranks used for classification and presented sorted by lineage (also available with `ganon report --sort lineage`):

```
unclassified                                                                  unclassified                                   0   0  0   2   2.02020
root           1        1                                                     root                                           0   0  97  97  97.97980
no rank        131567   1|131567                                              cellular organisms                             0   0  97  97  97.97980
superkingdom   2        1|131567|2                                            Bacteria                                       0   0  97  97  97.97980
phylum         1224     1|131567|2|1224                                       Proteobacteria                                 0   0  40  40  40.40404
class          1236     1|131567|2|1224|1236                                  Gammaproteobacteria                            0   0  12  12  12.12121
order          72274    1|131567|2|1224|1236|72274                            Pseudomonadales                                0   0  12  12  12.12121
family         468      1|131567|2|1224|1236|72274|468                        Moraxellaceae                                  0   0  12  12  12.12121
genus          469      1|131567|2|1224|1236|72274|468|469                    Acinetobacter                                  0   0  12  12  12.12121
species group  909768   1|131567|2|1224|1236|72274|468|469|909768             Acinetobacter calcoaceticus/baumannii complex  0   0  12  12  12.12121
species        470      1|131567|2|1224|1236|72274|468|469|909768|470         Acinetobacter baumannii                        12  0  0   12  12.12121
class          28211    1|131567|2|1224|28211                                 Alphaproteobacteria                            0   0  28  28  28.28283
order          204458   1|131567|2|1224|28211|204458                          Caulobacterales                                0   0  28  28  28.28283
family         76892    1|131567|2|1224|28211|204458|76892                    Caulobacteraceae                               0   0  28  28  28.28283
genus          75       1|131567|2|1224|28211|204458|76892|75                 Caulobacter                                    0   0  28  28  28.28283
species        366602   1|131567|2|1224|28211|204458|76892|75|366602          Caulobacter sp. K31                            28  0  0   28  28.28283
no rank        1783272  1|131567|2|1783272                                    Terrabacteria group                            0   0  57  57  57.57576
phylum         1239     1|131567|2|1783272|1239                               Firmicutes                                     0   0  57  57  57.57576
class          91061    1|131567|2|1783272|1239|91061                         Bacilli                                        0   0  57  57  57.57576
order          1385     1|131567|2|1783272|1239|91061|1385                    Bacillales                                     0   0  57  57  57.57576
family         186822   1|131567|2|1783272|1239|91061|1385|186822             Paenibacillaceae                               0   0  57  57  57.57576
genus          44249    1|131567|2|1783272|1239|91061|1385|186822|44249       Paenibacillus                                  0   0  57  57  57.57576
species        1406     1|131567|2|1783272|1239|91061|1385|186822|44249|1406  Paenibacillus polymyxa                         57  0  0   57  57.57576
```

## Hierarchical classification

Ganon classification can be performed in one or more databases at the same time. The databases can be provided in a hierarchical order. Multiple database classification can be performed providing several inputs for `--db-prefix`. They are required to be built with the same `--kmer-size` and `--window-size` values. To classify reads in a hierarchical order, `--hierarchy-labels` should be provided. `--abs-cutoff/--rel-cutoff` can be provided for each database. `--abs-filter/--rel-filter` can be provided for each hierarchy level. When using multiple hierarchical levels, output files will be generated for each level (use `--output-single` to generate a single output from multiple hierarchical levels).

For example: 

### Classifying reads against multiple databases:

```bash
ganon classify --db-prefix db1 db2 db3 \
               --rel-cutoff 0.75 \
               -r reads.fq.gz
```

Classification against 3 database (as if they were one) using the same error rate.

### Classifying reads against multiple databases with different error rates:

```bash
ganon classify --db-prefix db1 db2 db3 \
               --abs-cutoff  0   1   4 \
               -r reads.fq.gz
```

Classification against 3 database (as if they were one) using different error rates for each.

### Classifying reads against multiple databases hierarchically:

```bash
ganon classify --db-prefix            db1     db2      db3 \ 
               --hierarchy-labels 1_first 1_first 2_second \
               -r reads.fq.gz
```

In this example, reads are going to be classified first against db1 and db2. Reads without a valid match will be further classified against db3. `--hierarchy-labels` are strings and are going to be sorted to define the hierarchy order, disregarding input order.

### Classifying reads against multiple databases hierarchically with different error rates:

```bash
ganon classify --db-prefix            db1     db2      db3 \
               --hierarchy-labels 1_first 1_first 2_second \
               --rel-cutoff              1     0.5     0.25 \
               --abs-filter              0       1 \
               -r reads.fq.gz
```

In this example, classification will be performed with different error rates for each database. For each hierarchy (`1_first` and `2_second`) a different `--abs-filter` will be used.

## Choosing/explaining parameters

### ganon classify

#### --single-reads and/or --paired-reads

ganon accepts single-end or paired-end reads. In paired-end mode, reads are always reported with the header of the first pair. Paired-end reads are classified in a standard forward-reverse orientation. The max. read length accepted is 65535.

#### cutoff and filter (--abs-cutoff, --rel-cutoff, --abs-filter, --rel-filter)

ganon has two parameters to control a match between reads and references: `cutoff` and `filter`. Both can be used with absolute values (`--abs-cutoff` and `--abs-filter`) or relative values (`--rel-cutoff` and `--rel-filter`) interchangeably. They can also be disabled.

Every read can be classified against none, one or more references. What will be reported is the remaining matches after `cutoff` and `filter` thresholds are applied, based on the number of shared k-mers (or minimizers) between sequences.

The `cutoff` is the first. It should be set as a minimal value to consider a match between a read and a reference. Next the `filter` is applied to the remaining matches. `filter` thresholds are relative to the best scoring match. Here one can control how far from the best match we want to allow further matches. `cutoff` can be interpreted as the lower bound to discard spurious matches and `filter` as the fine tuning to what to keep.

For example:

Using `--kmer-size 19`, a certain read (100bp) has the following matches with the 5 references (`ref1..5`), sorted by shared k-mers:

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

since the `--rel-cutoff` threhsold is `82 * 0.25 = 21` (ceiling is applied). Further, with `--abs-filter 1`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.75 | --abs-filter 1 |
|-----------|---------------|-------------------|----------------|
| ref1      | 82            |                   |                |
| ref2      | 68            |                   |                |
| ~~ref3~~  | ~~44~~        |                   | X              |
| ~~ref4~~  | ~~25~~        |                   | X              |
| ~~ref5~~  | ~~20~~        | X                 |                |


since best match is 82 (0 errors), the filter parameter is allowing +1 error (the represents 63 to 81 shared k-mers). The error numbers are calculated using the q-gram lemma.

Currently, databases built with `--window-size` will only work with relative values. In this case, the relative values are not based on the maximum number of possible shared k-mers but on the actual number of minimizers extracted from the read.

Using absolute values for `cutoff` and `filter` is advised, since they will incorporate stricter limits defined by the q-gram lemma into your results. However, if your input consists of reads of different lengths, the use of `--rel-cutoff` is advised, since the absolute number of error in sequences of different sizes have different meanings. Note that in this case the use of `--abs-filter` is still possible and advised to better select matches based on the q-gram lemma.

A different `cutoff` can be set for every database in a multi or hierarchical database classification. A different `filter` can be set for every level of a hierarchical database classification.

Note that reads that remain with only one reference match (before or after `cutoff` and `filter` are applied) is considered an unique match.

#### minimizers or offset

`--window-size` or `--offset` can be used to speed-up analysis and reduce memory usage at the cost of sensitivity/precision on classification.

`--window-size` can be set with `ganon build` to activate the use of minimizers. It produces smaller database files and requires substantially less memory overall. It may increase building times but will have a huge benefit for classification times. We experience a 6x decrease in memory usage and 6x speed-up on classification (`--window-size 32 --kmer-size 19`). Sensitivity and precision were reduced by small margins. `--window-size` has to be greater or equal `--kmer-size`.

`--offset` can be set on (>1) and off (1) with `ganon classify` (for databases build without `--window-size`) . `--offset n` will only evaluate every nth k-mer of the input sequences. In our experiments, lower `--offset` values (e.g. 2,3,4) will drastically reduce running times on classification with very low impact on the sensitivity and precision of the results. `--offset` has to be smaller or equal `--kmer-size`.

### ganon build 

#### --max-filter-size and --bin-length 

The most useful variable to define the IBF size (.ibf file) is the `--max-filter-size`. It will set an approximate upper limit size for the file and estimate the `--bin-length` size based on it. However, there is a minimum size necessary to generate the filter given a set of references and chosen parameters. When using `--window-size`, the final filter size may be significantly smaller than calculated, depending on how similar are reference sequences in your dataset.

<details>
  <summary>More details</summary>

The IBF size is defined mainly on the amount of the input reference sequences (`-i`) but also can also be adjusted by a combination of parameters. Ganon will try to find the best `--bin-length` given `--max-fp`. Increasing `--max-fp` will generate smaller filters, but will generate more false positives in the classification step. If you know what you are doing, you can also directly set the size of the IBF with `--filter-size`.

`--bin-length` is the size in base pairs of each group for the taxonomic clustering (with TaxSBP). By default, `--fragment-length` will be the size of `--bin-length` - `--overlap-length`, meaning that sequences will be split with overlap to fit into the bins. For example: species X has 2 sequences of 120bp each. Considering `--bin-length 50` and `--overlap-length 10` (`--fragment-length 40` consequently) each of the sequences will be split into 50bp and put into a bin with overlap of 10bp, resulting in 3 bins for each sequence (6 in total for species X).

Such adjustment is necessary to equalize the size of each bin, since the IBF requires the individual bloom filters to be of the same size by definition. Building the IBF based on the biggest sequence group in your references will generate the lowest number of bins but a very sparse and gigantic IBF. Building the IBF based on the smallest sequence group in your references will generate the smallest IBF but with too many bins. A balance between those two is necessary to achieve small and fast filters.

</details>

#### --specialization

With the `--specialization` parameter is possible use an extra specialized "rank" as target for classification after taxonomic leaves. With this option, classification can be performed at strain, assembly, file, sequence level or any other custom level. `--specialization sequence` will use each each sequence (starting with ">") in the provided input files as a unique target to build the database. `--specialization file` will use each file as a target, useful for building at assembly/strain level if all sequences of the same group are in one file. `--specialization assembly` will retrieve assembly accessions from NCBI eutils and use them as target (this can take some time, since NCBI web services are limited on the amount of request per second). `--specialization custom` will use the 4th column of the `--seq-info-file` as target, allowing customized specializations.

### ganon update

#### --update-complete

By default `ganon update` will only add sequences provided with `--input-files` to an previously generated index. Using `--update-complete` it is possible to add and remove sequences from an index. When activating this option, ganon will consider that the files provided in `--input-files` are an actual representation of the index to build. It will automatically detect sequences that should be kept, inserted or removed given the input files and the information contained on the index to be updated.

#### --specialization

The same use of `--specialization` from `ganon build` (check above). `--specialization` can be changed in the update process (e.g. `--specialization assembly` was used to build and `--specialization file` is used to update). However, `--specialization` can only be used in `ganon update` if the database was built with it.

### ganon report

#### --report-type

By default, `ganon classify` and  `ganon report` generate a read-based report (`ganon report --report-type reads`) where each read classified is counted once, either to its unique or lca assignment. It is possible to generate the same report for all read matches (`ganon report --report-type matches`). In this case, multiple matches for each reads are reported to their targets (single or shared matches). Using `--report-type matches` will not show the unclassified number of reads and it will always sum up to 100% in the root node, since this reports the overall distribution of matches and not the amount of reads classified.

#### --split-hierarchy, --keep-hierarchy or --skip-hierarchy

When using multiple databases in different hierarchical levels to classify reads, it is possible to report them separately using `--split-hierarchy`. Once activated, one report will be generated for each hierarchical label. It is also possible to select or ignore specific hierarchical labels (e.g. for a label use for host removal) using `--keep-hierarchy` or `--skip-hierarchy`.

## Install without conda

<details>
  <summary>Instructions</summary>

### build dependencies

System packages:
- gcc >=7
- cmake >=3.10
- zlib

### run dependencies

System packages:
- python >=3.5
- pandas >=0.22.0
- gawk
- grep
- tar
- curl
- wget
- coreutils (zcat)

** Please make sure that the system packages are supported/installed in your environment. All other packages are installed in the next steps.

### Installing taxsbp, binpacking and pylca

```bash
git clone https://github.com/pirovc/pylca.git
cd pylca
git checkout d1474b2ec2c028963bafce278ccb69cc21c061fa #v1.0.0
python3 setup.py install
```

```bash
pip3 install binpacking==1.4.3
binpacking -h
```

```bash
git clone https://github.com/pirovc/taxsbp.git
cd taxsbp
git checkout 35ffb1e1a92f6199d757dfdd2f1971db29dd4070 # v1.1.1
python3 setup.py install
taxsbp -h
```

### Downloading and building ganon

```bash
git clone --recurse-submodules https://github.com/pirovc/ganon.git # ganon, catch2, cxxopts, seqan3
```
  
```bash
cd ganon
python3 setup.py install --record files.txt #optional
mkdir build_cpp
cd build_cpp
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF ..
make
sudo make install #optional
```

- to change install location (e.g. `/myprefix/bin/`), set the installation prefix in the cmake command with `-DCMAKE_INSTALL_PREFIX=/myprefix/ `

- use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.

If everything was properly installed, the following commands should show the help pages without errors:

```bash
ganon -h
```

### Run tests

```bash
python3 -m unittest discover -s tests/ganon/unit/
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/
cd build_cpp/
ctest -VV .
```

</details>

## Parameters

```
usage: ganon [-h] [-v] {build,update,classify,report,table} ...

ganon

positional arguments:
  {build,update,classify,report,table}
    build               Build ganon database
    update              Update ganon database
    classify            Classify reads
    report              Generate reports
    table               Generate table from reports

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
```

<details>
  <summary>ganon build</summary>

```
usage: ganon build [-h] -d DB_PREFIX [-i [INPUT_FILES ...]] [-t] [-r] [-m]
                   [-f] [-k] [-w] [--hash-functions] [--filter-size]
                   [--bin-length] [--fragment-length] [--overlap-length] [-s]
                   [--seq-info-mode [...]] [--seq-info-file]
                   [--taxdump-file [...]] [--input-directory]
                   [--input-extension] [--write-seq-info-file] [--verbose]
                   [--quiet]

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (.ibf, .map, .tax, .gnn will be
                        created)
  -i [INPUT_FILES ...], --input-files [INPUT_FILES ...]
                        Input reference sequence fasta files [.gz]

important arguments:
  -t , --threads        Number of sub-processes/threads to use. Default: 2
  -r , --rank           Target taxonomic rank for classification
                        [species,genus,...]. use "leaves" to use the leaf
                        taxonomic node assigned to each sequence as targets.
                        To use assembly, strain or further specializations,
                        check --specialization. Default: species
  -m , --max-filter-size 
                        Given an approx. upper limit in Megabytes (MB) for
                        filter/memory usage. When using --window-size, filter
                        may be significantly smaller after build depending on
                        the level of similarity of your input sequences.
                        [Mutually exclusive --bin-length]
  -f , --max-fp         Max. false positive rate for bloom filters [Mutually
                        exclusive --filter-size]. Default: 0.05

filter arguments:
  -k , --kmer-size      The k-mer size to split sequences. Default: 19
  -w , --window-size    The window-size to build filter with minimizers. 0 to
                        turn it off. Default: 0
  --hash-functions      The number of hash functions for the interleaved bloom
                        filter [0-5]. 0 to detect optimal value. Default: 0
  --filter-size         Fixed size for filter in Megabytes (MB) [Mutually
                        exclusive --max-fp]
  --bin-length          Maximum length (in bp) for each bin [Mutually
                        exclusive --max-filter-size]
  --fragment-length     Fragment sequences into specified length (in bp). Set
                        to 0 to not fragment sequences. Default: bin-length -
                        overlap-length
  --overlap-length      Fragment overlap length (in bp). Should be bigger than
                        the read length used for classification. Default: 300

other arguments:
  -s , --specialization 
                        Add extra specialized "rank" as target for
                        classification after taxonomic leaves. If set, --rank
                        is defaulted to leaves. Options:
                        [sequence,file,assembly,custom]. "sequence" will use
                        sequence accession as target. "file" uses the filename
                        as target. "assembly" will use assembly info from NCBI
                        as target. "custom" uses the 4th column of the file
                        provided in --seq-info-file as target.
  --seq-info-mode [ ...]
                        Automatic mode to retrieve tax. info and seq. length.
                        [auto,eutils] or one or more accession2taxid files
                        from NCBI [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot
                        dead_nucl dead_wgs dead_prot]. auto will either use
                        eutils for less than 50000 input sequences or nucl_gb
                        nucl_wgs. Alternatively a file can be directly
                        provided (see --seq-info-file). Default: auto
  --seq-info-file       Tab-separated file with sequence information (seqid
                        <tab> seq.len <tab> taxid [<tab> specialization])
                        [Mutually exclusive --seq-info-mode].
  --taxdump-file [ ...]
                        Force use of a specific version of the
                        (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
                        file(s) from NCBI Taxonomy (otherwise it will be
                        automatically downloaded).
  --input-directory     Directory containing input files
  --input-extension     Extension of files to use with --input-directory
                        (provide it without * expansion, e.g. ".fna.gz")
  --write-seq-info-file
                        Write sequence information to DB_PREFIX.seqinfo.txt
  --verbose             Verbose output mode
  --quiet               Quiet output mode
```

</details>

<details>
  <summary>ganon update</summary>

```
usage: ganon update [-h] -d DB_PREFIX [-i [INPUT_FILES ...]] [-o] [-t] [-s]
                    [--seq-info-mode [...]] [--seq-info-file]
                    [--taxdump-file [...]] [--input-directory]
                    [--input-extension] [--update-complete]
                    [--write-seq-info-file] [--verbose] [--quiet]

required arguments:
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database input prefix (.ibf, .map, .tax, .gnn)
  -i [INPUT_FILES ...], --input-files [INPUT_FILES ...]
                        Input reference sequence fasta files [.gz] to be
                        included to the database. Complete set of updated
                        sequences should be provided when using --update-
                        complete

optional arguments:
  -h, --help            show this help message and exit
  -o , --output-db-prefix 
                        Output database prefix (.ibf, .map, .tax, .gnn).
                        Default: overwrite current --db-prefix
  -t , --threads        Number of sub-processes/threads to use. Default: 2
  -s , --specialization 
                        Change specialization mode. Can only be used if
                        database was built with some specialization. Options:
                        [sequence,file,assembly,custom]. "sequence" will use
                        sequence accession as target. "file" uses the filename
                        as target. "assembly" will use assembly info from NCBI
                        as target. "custom" uses the 4th column of the file
                        provided in --seq-info-file as target.
  --seq-info-mode [ ...]
                        Automatic mode to retrieve tax. info and seq. length.
                        [auto,eutils] or one or more accession2taxid files
                        from NCBI [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot
                        dead_nucl dead_wgs dead_prot]. auto will either use
                        eutils for less than 50000 input sequences or nucl_gb
                        nucl_wgs. Alternatively a file can be directly
                        provided (see --seq-info-file). Default: auto
  --seq-info-file       Tab-separated file with sequence information (seqid
                        <tab> seq.len <tab> taxid [<tab> assembly id])
                        [Mutually exclusive --seq-info]
  --taxdump-file [ ...]
                        Force use of a specific version of the
                        (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
                        file(s) from NCBI Taxonomy (otherwise it will be
                        automatically downloaded)
  --input-directory     Directory containing input files
  --input-extension     Extension of files to use with --input-directory
                        (provide it without * expansion, e.g. ".fna.gz")
  --update-complete     Update adding and removing sequences. Input files
                        should represent the complete updated set of
                        references, not only new sequences.
  --write-seq-info-file
                        Write sequence information to DB_PREFIX.seqinfo.txt
  --verbose             Verbose output mode
  --quiet               Quiet output mode
```

</details>

<details>
  <summary>ganon classify</summary>

```
usage: ganon classify [-h] -d [DB_PREFIX ...] [-s [reads.fq[.gz] ...]]
                      [-p [reads.1.fq[.gz] reads.2.fq[.gz] ...]] [-c [...]]
                      [-b [...]] [-e [...]] [-a [...]] [-o] [--output-lca]
                      [--output-all] [--output-unclassified] [--output-single]
                      [-l [...]] [-f] [-t] [-r [...]] [--verbose] [--quiet]

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -d [DB_PREFIX ...], --db-prefix [DB_PREFIX ...]
                        Database input prefix[es]
  -s [reads.fq[.gz] ...], --single-reads [reads.fq[.gz] ...]
                        Multi-fastq[.gz] file[s] to classify
  -p [reads.1.fq[.gz] reads.2.fq[.gz] ...], --paired-reads [reads.1.fq[.gz] reads.2.fq[.gz] ...]
                        Multi-fastq[.gz] pairs of file[s] to classify

cutoff/filter arguments:
  -c [ ...], --rel-cutoff [ ...]
                        Min. relative percentage of k-mers necessary to
                        consider a match. Generally used to cutoff low
                        similarity matches. Single value or one per database
                        (e.g. 0.5 0.7 1 0.25). 0 for no cutoff. [Mutually
                        exclusive --rel-cutoff] Default: 0.25
  -b [ ...], --abs-cutoff [ ...]
                        Max. absolute number of errors allowed to consider a
                        match. Generally used to cutoff low similarity
                        matches. Single value or one per database (e.g. 3 3 4
                        0). -1 for no cutoff. [Mutually exclusive --abs-
                        cutoff]
  -e [ ...], --rel-filter [ ...]
                        Additional relative percentage of k-mers (relative to
                        the best match) to keep a match (applied after
                        cutoff). Single value or one per hierarchy (e.g. 0.1 0
                        0.25). 1 for no filter. [Mutually exclusive --abs-
                        filter]
  -a [ ...], --abs-filter [ ...]
                        Additional absolute number of errors (relative to the
                        best match) to keep a match (applied after cutoff).
                        Single value or one per hierarchy (e.g. 0 2 1). -1 for
                        no filter. [Mutually exclusive --rel-filter] Default:
                        0

output arguments:
  -o , --output-prefix 
                        Output prefix to print report (.rep). Empty to output
                        to STDOUT
  --output-lca          Output an additional file with one lca match for each
                        read (.lca)
  --output-all          Output an additional file with all matches. File can
                        be very large (.all)
  --output-unclassified
                        Output an additional file with unclassified read
                        headers (.unc)
  --output-single       When using multiple hierarchical levels, output
                        everything in one file instead of one per hierarchy

other arguments:
  -l [ ...], --hierarchy-labels [ ...]
                        Hierarchy definition, one for each database input. Can
                        also be a string, but input will be sorted to define
                        order (e.g. 1 1 2 3). Default: H1
  -f , --offset         Number of k-mers to skip during classification. Can
                        speed up analysis but may reduce recall. (e.g. 1 = all
                        k-mers, 3 = every 3rd k-mer). Default: 1
  -t , --threads        Number of sub-processes/threads to use. Default: 3
  -r [ ...], --ranks [ ...]
                        Ranks to show in the report (.tre). "all" for all
                        identified ranks. empty for default ranks:
                        superkingdom phylum class order family genus species
                        assembly. This file can be re-generated with the ganon
                        report command.
  --verbose             Verbose output mode
  --quiet               Quiet output mode
```

</details>

<details>
  <summary>ganon report</summary>

```
usage: ganon report [-h] [--min-count] [--max-count] [--names [...]]
                    [--names-with [...]] [--taxids [...]] [-i [REP_FILES ...]]
                    -o OUTPUT_PREFIX [-d [...]] [-f] [-e] [-r [...]] [-s] [-y]
                    [-p [...]] [-k [...]] [--taxdump-file [...]]
                    [--input-directory] [--input-extension] [--verbose]
                    [--quiet]

required arguments:
  -i [REP_FILES ...], --rep-files [REP_FILES ...]
                        One or more *.rep files from ganon classify
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix for report file "{output_prefix}.tre".
                        In case of multiple files, the base input filename
                        will be appended at the end of the output file
                        "{output_prefix + FILENAME}.tre"

filter arguments:
  --min-count           Minimum number/percentage of counts to keep an taxa
                        [values between 0-1 for percentage, >1 specific
                        number]
  --max-count           Maximum number/percentage of counts to keep an taxa
                        [values between 0-1 for percentage, >1 specific
                        number]
  --names [ ...]        Show only entries matching exact names of the provided
                        list
  --names-with [ ...]   Show entries containing full or partial names of the
                        provided list
  --taxids [ ...]       One or more taxids to report (including children taxa)

optional arguments:
  -h, --help            show this help message and exit
  -d [ ...], --db-prefix [ ...]
                        Database prefix[es] used for classification (in any
                        order). Only ".tax" file is required. If not provided,
                        new taxonomy will be downloaded
  -f , --output-format 
                        Output format [text, tsv, csv]. text outputs a
                        tabulated formatted text file for better
                        visualization. Default: tsv
  -e , --report-type    Type of report to generate [reads, matches]. Default:
                        reads
  -r [ ...], --ranks [ ...]
                        Ranks to report ["", "all", custom list] "all" for all
                        possible ranks. empty for default ranks (superkingdom
                        phylum class order family genus species assembly).
                        Default: ""
  -s , --sort           Sort report by [rank, lineage, count, unique].
                        Default: rank (with custom --ranks) or lineage (with
                        --ranks all)
  -y, --split-hierarchy
                        Split output reports by hierarchy (from ganon classify
                        --hierarchy-labels). If activated, the output files
                        will be named as "{output_prefix}.{hierarchy}.tre"
  -p [ ...], --skip-hierarchy [ ...]
                        One or more hierarchies to skip in the report (from
                        ganon classify --hierarchy-labels)
  -k [ ...], --keep-hierarchy [ ...]
                        One or more hierarchies to keep in the report (from
                        ganon classify --hierarchy-labels)
  --taxdump-file [ ...]
                        Force use of a specific version of the
                        (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
                        file(s) from NCBI Taxonomy (otherwise it will be
                        automatically downloaded)
  --input-directory     Directory containing input files
  --input-extension     Extension of files to use with --input-directory
                        (provide it without * expansion, e.g. ".rep")
  --verbose             Verbose output mode
  --quiet               Quiet output mode
```

</details>

<details>
  <summary>ganon table</summary>

```
usage: ganon table [-h] [--min-count] [--max-count] [--names [...]]
                   [--names-with [...]] [--taxids [...]] [-i [TRE_FILES ...]]
                   -o OUTPUT_FILE [-l] [-f] [-t] [-a] [-m] [-r] [-n]
                   [--header] [--unclassified-label] [--filtered-label]
                   [--skip-zeros] [--transpose] [--input-directory]
                   [--input-extension] [--verbose] [--quiet]

required arguments:
  -i [TRE_FILES ...], --tre-files [TRE_FILES ...]
                        Report files (.tre) from ganon classify/report to make
                        the table
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output filename for the table

filter arguments:
  --min-count           Minimum number/percentage of counts to keep an taxa
                        [values between 0-1 for percentage, >1 specific
                        number]
  --max-count           Maximum number/percentage of counts to keep an taxa
                        [values between 0-1 for percentage, >1 specific
                        number]
  --names [ ...]        Show only entries matching exact names of the provided
                        list
  --names-with [ ...]   Show entries containing full or partial names of the
                        provided list
  --taxids [ ...]       One or more taxids to report (including children taxa)

optional arguments:
  -h, --help            show this help message and exit
  -l , --output-value   Output value on the table [percentage, counts].
                        percentage values are reported between [0-1]. Default:
                        counts
  -f , --output-format 
                        Output format [tsv, csv]. Default: tsv
  -t , --top-sample     Top hits of each sample individually
  -a , --top-all        Top hits of all samples (ranked by percentage)
  -m , --min-frequency 
                        Minimum number/percentage of files containing an taxa
                        to keep the taxa [values between 0-1 for percentage,
                        >1 specific number]
  -r , --rank           Define specific rank to report. Empty will report all
                        ranks.
  -n, --no-root         Do not report root node entry and lineage. Direct and
                        shared matches to root will be accounted as
                        unclassified
  --header              Header information [name, taxid, lineage]. Default:
                        name
  --unclassified-label 
                        Add column with unclassified count/percentage with the
                        chosen label. May be the same as --filtered-label
                        (e.g. unassigned)
  --filtered-label      Add column with filtered count/percentage with the
                        chosen label. May be the same as --unclassified-label
                        (e.g. unassigned)
  --skip-zeros          Do not print lines with only zero count/percentage
  --transpose           Transpose output table (taxa as cols and files as
                        rows)
  --input-directory     Directory containing input files
  --input-extension     Extension of files to use with --input-directory
                        (provide it without * expansion, e.g. ".tre")
  --verbose             Verbose output mode
  --quiet               Quiet output mode
```

</details>
