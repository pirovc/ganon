# Output files

## ganon build/build-custom/update

Every run on `ganon build`, `ganon build-custom` or `ganon update` will generate the following database files:

 - {prefix}**.ibf/.hibf**: main bloom filter index file, extension based on the `--filter-type` option.
 - {prefix}**.tax**: taxonomy tree, only generated if `--taxonomy` is used *(fields: target/node, parent, rank, name, genome size)*.
 - {prefix}**_files/**: (`ganon build` only) folder containing downloaded reference sequence and auxiliary files. Not necessary for classification. Keep this folder if the database will be update later. Otherwise it can be deleted.

!!! warning
    Database files generated with version 1.2.0 or higher are not compatible with older versions.

## ganon classify
 
- {prefix}**.tre**: full report file (see below)
- {prefix}**.rep**: plain report of the run with only targets that received a match. Can be used to re-generate full reports (.tre) with `ganon report`. At the end prints 2 extra lines with `#total_classified` and `#total_unclassified`. Fields
    - 1: hierarchy label
    - 2: target
    - 3: # total matches
    - 4: # unique reads
    - 5: # lca reads
    - 6: rank
    - 7: name
- {prefix}**.one**: output with one match for each classified read after EM or LCA algorithm. Only generated with `--output-one` active. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierarchy}.one *(fields: read identifier, target, (max) k-mer/minimizer count)*
- {prefix}**.all**: output with all matches for each read. Only generated with `--output-all` active **Warning: file can be very large**. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierarchy}.all *(fields: read identifier, target, k-mer/minimizer count)*

## ganon report

- {prefix}**.tre**: tab-separated tree-like report with cumulative counts and taxonomic lineage. There are several possible `--report-type`. More information on the different types of reports can be found [here](#report-type---report-type):
    - *abundance*: will attempt to estimate **taxonomic abundances** by re-disributing read counts from LCA matches and correcting sequence abundance by approximate genome sizes.
    - *reads*: **sequence abundances**, reports the proportion of sequences assigned to a taxa, each read classified is counted once.
    - *dist*: like *reads* with read count re-distribution.
    - *corr*: like *reads* with correction by genome size.
    - *matches*: every match is reported to their original target, including multiple and shared matches.

Each line in this report is a taxonomic entry (including the root node), with the following fields: 

| col | field        | obs                                                                                                                                                                                                                            | example                                       |
|-----|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------|
| 1   | rank         |                                                                                                                                                                                                                                | phylum                                        |
| 2   | target       | taxonomic id. or specialization (assembly id.)                                                                                                                                                                                 | 562                                           |
| 3   | lineage      |                                                                                                                                                                                                                                | 1\|131567\|2\|1224\|28211\|766\|942\|768\|769 |
| 4   | name         |                                                                                                                                                                                                                                | Chromobacterium rhizoryzae                    |
| 5   | # unique     | number of reads that matched exclusively to this target                                                                                                                                                                        | 5                                             |
| 6   | # shared     | number of reads with non-unique matches directly assigned to this target. Represents the LCA matches (`--report-type reads`), re-assigned matches (`--report-type abundance/dist`) or shared matches (`--report-type matches`) | 10                                            |
| 7   | # children   | number of unique and shared assignments to all children nodes of this target                                                                                                                                                   | 20                                            |
| 8   | # cumulative | the sum of the unique, shared and children assignments up-to this target                                                                                                                                                       | 35                                            |
| 9   | % cumulative | percentage of assignments or estimated relative abundance for `--report-type abundance`                                                                                                                                                                     | 43.24                                         |

- The first line of the report file will show the number of unclassified reads (not for `--report-type matches` or `--normalize`)

- The CAMI challenge [bioboxes profiling format](https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd) is supported using `--output-format bioboxes`. In this format, only values for the percentage/abundance (col. 9) are reported. The root node and unclassified entries are omitted.

- The sum of cumulative assignments for the unclassified and root lines is 100%. The final cumulative sum of reads/matches may be under 100% if any filter is successfully applied and/or hierarchical selection is selected (keep/skip/split).

- For all report type but `matches`, only taxa that received direct read matches, either unique or by LCA assignment, are considered. Some reads may have only shared matches and will not be reported directly but will be accounted for on some parent level. To visualize those matches, create a report with `--report-type matches` or use directly the file {prefix}**.rep**.

## ganon table

 - {output_file}: a tab-separated file with counts/percentages of taxa for multiple samples
 
---

<details>
  <summary>Examples of output files</summary>

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

with `--output-format bioboxes` 

```
@Version:0.10.0
@SampleID:example.rep H1
@Ranks:superkingdom|phylum|class|order|family|genus|species|assembly
@Taxonomy:db.tax
@@TAXID  RANK          TAXPATH   TAXPATHSN                 PERCENTAGE
2        superkingdom  2         Bacteria                  100.00000
1224     phylum        2|1224    Bacteria|Proteobacteria   56.89782
201174   phylum        2|201174  Bacteria|Actinobacteria   21.84869
1239     phylum        2|1239    Bacteria|Firmicutes       9.75197
976      phylum        2|976     Bacteria|Bacteroidota     6.15297
1117     phylum        2|1117    Bacteria|Cyanobacteria    2.23146
203682   phylum        2|203682  Bacteria|Planctomycetota  1.23353
57723    phylum        2|57723   Bacteria|Acidobacteria    0.52549
200795   phylum        2|200795  Bacteria|Chloroflexi      0.31118
```
</details>
<br>