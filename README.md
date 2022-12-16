# ganon [![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html) [![Publication](https://img.shields.io/badge/DOI-10.1101%2F406017-blue)](https://dx.doi.org/10.1093/bioinformatics/btaa458)

ganon classifies short DNA sequences against large sets of genomic reference sequences efficiently. It automatically downloads, builds and updates commonly used repositories (refseq/genbank), performs taxonomic (ncbi or gtdb) and hierarchical classification, generates taxonomic and/or sequence abundance reports, generates contingency tables among many other [features](#Features).

---

## Quick install/usage guide

### Install with conda

```sh
conda install -c bioconda -c conda-forge ganon
```

### Download and build
```bash
# Archaeal complete genome sequences from NCBI RefSeq
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 12
```

### Classify
```bash
ganon classify --db-prefix arc_cg_rs --output-prefix classify_results --single-reads my_reads.fq.gz --threads 12
```

### Re-generate reports and create tables from multiple reports
```bash
ganon report --db-prefix arc_cg_rs --input classify_results.rep --output-prefix filtered_report --min-count 0.01
ganon table --input classify_results.tre filtered_report.tre --output-file output_table.tsv --top-sample 10
```

### Update the database at a later time point
```bash
ganon update --db-prefix arc_cg_rs --threads 12
```

---

- [Details](#details)
- [Examples](#examples)
- [Output files](#output-files)
- [Building customized databases](#building-customized-databases)
- [Multiple and Hierarchical classification](#multiple-and-hierarchical-classification)
- [Choosing and explaining parameters](#choosing-and-explaining-parameters)
- [Parameters](#parameters)

## Details

ganon is designed to index large sets of genomic reference sequences and to classify short reads against them efficiently. The tool uses Interleaved Bloom Filters as indices based on k-mers/minimizers. It was mainly developed, but not limited, to the metagenomics classification problem: quickly assign short fragments to their closest reference among thousands of references. After classification, taxonomic abundance is estimated and reported.

### Features

- NCBI and GTDB native support for taxonomic classification (but also runs without taxonomy)
- integrated download of commonly used reference sequences from RefSeq/Genbank (`ganon build`)
- update indices incrementally (`ganon update`)
- customizable build for pre-downloaded or non-standard sequence files (`ganon build-custom`)
- build and classify at different taxonomic levels, file, sequence, strain/assembly or custom specialization
- perform [hierarchical classification](#multiple-and-hierarchical-classification): use several databases in any order
- [report](#report) the lowest common ancestor (LCA) but also multiple and unique matches for every read
- [report](#report) sequence or taxonomic abundances as well as total number of matches
- generate reports and contingency tables for multi-sample studies with several filter options

ganon achieved very good results in [our own evaluations](https://dx.doi.org/10.1093/bioinformatics/btaa458) but also in independent evaluations: [LEMMI](https://lemmi-v1.ezlab.org/), [LEMMI v2](https://lemmi.ezlab.org/) and [CAMI2](https://dx.doi.org/10.1038/s41592-022-01431-4)

### Installation guide

The easiest way to install ganon is via conda, using the bioconda and conda-forge channels:

```bash
conda install -c bioconda -c conda-forge ganon
```

However, there are possible performance benefits compiling ganon from source in the target machine rather than using the conda version. To do so, please follow the instructions below:

<details>
  <summary>Instructions</summary>

#### build dependencies

System packages:
- gcc >=7
- cmake >=3.10
- zlib

#### run dependencies

System packages:
- python >=3.6
- pandas >=1.1.0
- [multitax](https://github.com/pirovc/multitax) >=1.2.1

```bash
python3 -V # >=3.6
python3 -m pip install "pandas>=1.1.0" "multitax>=1.2.1"
```

#### Downloading and building ganon + submodules

```bash
git clone --recurse-submodules https://github.com/pirovc/ganon.git
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

#### Run tests

```bash
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/ #optional - downloads large files
cd build_cpp/
ctest -VV .
```
</details>

## Examples

### Commonly used reference set for metagenomics analysis (archaea+bacteria+fungi+viral from NCBI RefSeq, complete genomes, one assembly per taxa)

```bash
ganon build --db-prefix abfv_rs_cg_t1a --organism-group archaea bacteria fungi viral --source refseq --taxonomy ncbi --complete-genomes --threads 12 --top 1
```

### GTDB database, one assembly per taxa
```bash
ganon build --db-prefix complete_gtdb --organism-group archaea bacteria --source refseq genbank --taxonomy gtdb --threads 12 --top 1
```

### Database based on specific taxonomic identifiers (203492 - Fusobacteriaceae)
```bash
# NCBI
ganon build --db-prefix fuso_ncbi --taxid "203492" --source refseq genbank --taxonomy ncbi --threads 12
# GTDB
ganon build --db-prefix fuso_gtdb --taxid "f__Fusobacteriaceae" --source refseq genbank --taxonomy gtdb --threads 12
```

### Customized database at assembly level 
```bash
ganon build-custom --db-prefix my_db --input my_big_fasta_file.fasta.gz --level assembly --threads 12
```

### Customized database at species level build based on files previously downloaded with [genome_updater](https://github.com/pirovc/genome_updater)
```bash
ganon build-custom --db-prefix custom_db_gu --input myfiles/2022-06-28_10-02-14/files/ --level species --ncbi-file-info outfolder/2022-06-28_10-02-14/assembly_summary.txt --threads 12
```

### Customized database with sequence as target (to classify reads at sequence level)
```bash
ganon build-custom --db-prefix seq_target --input myfiles/2022-06-28_10-02-14/files/ --ncbi-file-info outfolder/2022-06-28_10-02-14/assembly_summary.txt --input-target sequence --threads 12
```
## Output files

### build/update

Every run on `ganon build`, `ganon build-custom` or `ganon update` will generate the following database files:

 - {prefix}**.ibf**: main interleaved bloom filter index file
 - {prefix}**.tax**: taxonomy tree, only generated if `--taxonomy` is used *(fields: target/node, parent, rank, name, genome size)*.
 - {prefix}**_files/**: (`ganon build` only) folder containing downloaded reference sequence and auxiliary files. Not necessary for classification. Keep this folder if the database will be update later. Otherwise it can be deleted.

*Obs: Database files generated with version 1.2.0 or higher are not compatible with older versions.*

### classify
 
 - {prefix}**.tre**: full report file (see below)
 - {prefix}**.rep**: plain report of the run with only targets that received a match. Can be used to re-generate full reports (.tre) with `ganon report`. At the end prints 2 extra lines with `#total_classified` and `#total_unclassified`. Fields:
   - 1: hierarchy label
   - 2: target
   - 3: # total matches
   - 4: # unique reads
   - 5: # lca reads
   - 6: rank
   - 7: name
 - {prefix}**.lca**: output with one match for each classified read after LCA. Only generated with `--output-lca` active. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierarchy}.lca *(fields: read identifier, target, (max) k-mer/minimizer count)*
 - {prefix}**.all**: output with all matches for each read. Only generated with `--output-all` active **Warning: file can be very large**. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierarchy}.all *(fields: read identifier, target, k-mer/minimizer count)*

### report

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

- The first line of the report file will show the number of unclassified reads (not for `--report-type matches`)

- The CAMI challenge [bioboxes profiling format](https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd) is supported using `--output-format bioboxes`. In this format, only values for the percentage/abundance (col. 9) aer reported. The root node and unclassified entries are ommited.

- The sum of cumulative assignments for the unclassified and root lines is 100%. The final cumulative sum of reads/matches may be under 100% if any filter is successfully applied and/or hierarchical selection is selected (keep/skip/split).

- For all report type but `matches`, only taxa that received direct read matches, either unique or by LCA assignment, are considered. Some reads may have only shared matches and will not be reported directly but will be accounted for on some parent level. To visualize those matches, create a report with `--report-type matches` or use directly the file {prefix}**.rep**.

### table

 - {output_file}: a tab-separated file with counts/percentages of taxa for multiple samples
 
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

## Building customized databases

Besides the automated download and build (`ganon build`) ganon provides a highly customizable build procedure (`ganon build-custom`) to create databases. 

To use custom sequences, just provide them with `--input`. ganon will try to retrieve all necessary information necessary to build a database.

*ganon expects assembly accessions if building by file (e.g. filename should be similar as `GCA_002211645.1_ASM221164v1_genomic.fna.gz`) or accession version if building by sequence (e.g. headers should look like `>CP022124.1 Fusobacterium nu...`).* More information about building by file or sequence can be found [here](#target-file-or-sequence---input-target).

It is also possible to use non-standard accessions and headers to build databases with `--input-file`. This file should contain the following fields (tab-separated): file, [target, node, specialization, specialization_name].

<details>
  <summary>Examples of --input-file</summary>

Using `--input-target sequence`:

```
sequences.fasta HEADER1
sequences.fasta HEADER2
sequences.fasta HEADER3
others.fasta HEADER4
others.fasta HEADER5
```

or using `--input-target file`:

```
sequences.fasta FILE_A
others.fasta FILE_B
```

Nodes can be provided to link the data with taxonomy. For example (using `--taxonomy ncbi`):

```
sequences.fasta HEADER1 562
sequences.fasta HEADER2 562
sequences.fasta HEADER3 562
others.fasta HEADER4  623
others.fasta HEADER5  623
```

Specialization can be used to create a additional classification level after the taxonomic leaves. For example (using `--level custom`):

```
sequences.fasta HEADER1 562 ID443 Escherichia coli TW10119
sequences.fasta HEADER2 562 ID297 Escherichia coli PCN079
sequences.fasta HEADER3 562 ID8873  Escherichia coli P0301867.7
others.fasta HEADER4  623 ID2241  Shigella flexneri 1a
others.fasta HEADER5  623 ID4422  Shigella flexneri 1b
```
</details>

## Multiple and Hierarchical classification

`ganon classify` can be performed in multiple databases at the same time. The databases can also be provided in a hierarchical order. 

Multiple database classification can be performed providing several inputs for `--db-prefix`. They are required to be built with the same `--kmer-size` and `--window-size` values. Multiple databases are considered as one (as if built together) and redundancy in content (same reference in two or more databases) is allowed.

To classify reads in a hierarchical order, `--hierarchy-labels` should be provided. When using multiple hierarchical levels, output files will be generated for each level (use `--output-single` to generate a single output from multiple hierarchical levels). Please note that some parameters are set for each database (e.g. `--rel-cutoff`) while others are set for each hierarchical level (e.g. `--rel-filter`)

<details>
  <summary>Examples</summary>

### Classifying reads against multiple databases:

```bash
ganon classify --db-prefix db1 db2 db3 \
               --rel-cutoff 0.75 \
               --single-reads reads.fq.gz
```

Classification against 3 database (as if they were one) using the same cutoff.

### Classifying reads against multiple databases with different cutoffs:

```bash
ganon classify --db-prefix  db1 db2 db3 \
               --rel-cutoff 0.2 0.3 0.1 \
               --single-reads reads.fq.gz
```

Classification against 3 database (as if they were one) using different error rates for each.

### Classifying reads against multiple databases hierarchically:

```bash
ganon classify --db-prefix            db1     db2      db3 \
               --hierarchy-labels 1_first 1_first 2_second \
               --single-reads reads.fq.gz
```

In this example, reads are going to be classified first against db1 and db2. Reads without a valid match will be further classified against db3. `--hierarchy-labels` are strings and are going to be sorted to define the hierarchy order, disregarding input order.

### Classifying reads against multiple databases hierarchically with different cutoffs:

```bash
ganon classify --db-prefix            db1     db2      db3 \
               --hierarchy-labels 1_first 1_first 2_second \
               --rel-cutoff             1     0.5     0.25 \
               --rel-filter           0.1              0.5 \
               --single-reads reads.fq.gz
```

In this example, classification will be performed with different `--rel-cutoff` for each database. For each hierarchy levels (`1_first` and `2_second`) a different `--rel-filter` will be used.

</details>

## Choosing and explaining parameters

The most important parameters and trade-offs are:

- `ganon build` `--window-size --kmer-size`: the *window* value should always be same or larger than *kmer* value. The larger the difference between them, the smaller the database size will be. However, some loss of sensitivity/precision is expected when using the database for classification. Larger kmer values (e.g. `31`) will improve classification, specially read binning, at a cost of larger database sizes.
- `ganon classify` `--rel-cutoff`: will define the matches between reads and database. Higher `--rel-cutoff` values will improve precision with some loss of sensitivity with a decrease in unique matches but an increase in overall matches. For taxonomic profiling, a higher value between `0.4` and `0.8` may provide better results. For using the read classification/binning output, lower values between `0.2` and `0.4` are recommended.
- `ganon classify` `--rel-filter`: filter matches after cutoff is applied. Usually set between `0` and `0.2`. With higher filter values, further good scoring matches will be reported.
- `ganon report` `--report-type`: reports either taxonomic, sequence or matches abundances. Use `corr` or `abundance` for taxonomic profiling, `reads` or `dist` for sequence profiling and `matches` to report a summary of all matches.
- `ganon report` `--min-count`: cutoff to discard underrepresented taxa. Useful to remove the common long tail of spurious matches and false positives when performing classification. Values between `0.0001` (0.01%) and `0.001` (0.1%) improved sensitivity and precision in our evaluations. The higher the value, the more precise the outcome, with a sensitivity loss. Alternatively `--top-percentile` can be used to keep a relative amount of taxa instead a hard cutoff.

The numeric values above are averages from several experiments with different sample types and database contents. They may not work as expected for your data. If you are not sure which values to use or see something unexpected, please open an [issue](https://github.com/pirovc/ganon/issues).

Below, a more in-depth explanation of important parameters:

### ganon build 

#### filter false positive and size (--max-fp, --filter-size)

ganon indices are based on bloom filters and can have false positive matches. This can be controlled with `--max-fp` parameter. The lower the `--max-fp`, the less chances of false positives, but the larger the database size will be. For example, with `--max-fp 0.01` the database will be build so any target (e.g. assembly, specificed with `--level`) will have 1 in a 100 change of reporting a false match (between minimizer k-mers).

Alternatively, one can set a specific size for the final index with `--filer-size`. When using this option, please observe the theoretic false positive of the index reported at the end of the building process.

#### minimizers (--window-size, --kmer-size)

in `ganon build`, when `--window-size` > `--kmer-size` minimizers are used. That means that for a every window, a single k-mer will be selected. It produces smaller database files and requires substantially less memory overall. It may increase building times but will have a huge benefit for classification times. Sensitivity and precision can be reduced by small margins. If `--window-size` = `--kmer-size`, all k-mers are going to be used to build the database.

### ganon classify

#### reads (--single-reads, --paired-reads)

ganon accepts single-end and paired-end reads. Both types can be use at the same time. In paired-end mode, reads are always reported with the header of the first pair. Paired-end reads are classified in a standard forward-reverse orientation. The max. read length accepted is 65535 (accounting both reads in paired mode).

#### cutoff and filter (--rel-cutoff, --rel-filter)

ganon has two parameters to control a match between reads and references: `--rel-cutoff` and `--rel-filter`.

Every read can be classified against none, one or more references. What will be reported is the remaining matches after `cutoff` and `filter` thresholds are applied, based on the number of shared minimizers (or k-mers) between sequences.

The `cutoff` is the first. It should be set as a minimal value to consider a match between a read and a reference. Next the `filter` is applied to the remaining matches. `filter` thresholds are relative to the best scoring match and control how far from the best match further matches are allowed. `cutoff` can be interpreted as the lower bound to discard spurious matches and `filter` as the fine tuning to control what to keep.

<details>
  <summary>Example</summary>

Using `--kmer-size 19` (and `--window-size 19` to simplify the example), a certain read (100bp) has the following matches with the 5 references (`ref1..5`), sorted by shared k-mers:

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

</details>

For databases built with `--window-size`, the relative values are not based on the maximum number of possible shared k-mers but on the actual number of unique minimizers extracted from the read.

A different `cutoff` can be set for every database in a multiple or hierarchical database classification. A different `filter` can be set for every level of a hierarchical database classification.

Note that reads that remain with only one reference match (after `cutoff` and `filter` are applied) are considered a unique match.

### ganon build-custom

#### Target file or sequence (--input-target) 

Customized builds can be done either by file or sequence. `--input-target file` will consider every file provided with `--input` a single unit. `--input-target sequence` will use every sequence as a unit.

`--input-target file` is the default behavior and most efficient way to build databases. `--input-target sequence` should only be used when the input sequences are stored in a single file or when classification at sequence level is desired.

#### Build level (--level)

The `--level` parameter defines the max. depth of the database for classification. This parameter is relevant because the `--max-fp` is going to be guaranteed at the `--level` chosen. By default, the level will be the same as `--input-target`, meaning that classification will be done either at file or sequence level.

Alternatively, `--level assembly` will link the file or sequence target information with assembly accessions retrieved from NCBI servers. `--level leaves` or `--level species` (or genus, family, ...) will link the targets with taxonomic information and prune the tree at the chosen level. `--level custom` will use specialization level define in the `--input-file`.

#### Genome sizes (--genome-size-files)

Ganon will automatically download auxiliary files to define an approximate genome size for each entry in the taxonomic tree. For `--taxonomy ncbi` the [species_genome_size.txt.gz](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/) is used. For `--taxonomy gtdb` the [\*_metadata.tar.gz](https://data.gtdb.ecogenomic.org/releases/latest/) files are used. Those files can be directly provided with the `--genome-size-files` argument.

Genome sizes of parent nodes are calculated as the average of the respective children nodes. Other nodes without direct assigned genome sizes will use the closest parent with a pre-calculated genome size. The genome sizes are stored in the [ganon database](#buildupdate).

#### Retrieving info (--ncbi-sequence-info, --ncbi-file-info)

Further taxonomy and assembly linking information has to be collected to properly build the database. `--ncbi-sequence-info` and `--ncbi-file-info` allow customizations on this step.

When `--input-target sequence`, `--ncbi-sequence-info` argument allows the use of NCBI e-utils webservices (`eutils`) or downloads accession2taxid files to extract target information (options `nucl_gb` `nucl_wgs` `nucl_est` `nucl_gss` `pdb` `prot` `dead_nucl` `dead_wgs` `dead_prot`). By default, ganon uses `eutils` up-to 50000 input sequences, otherwise it downloads `nucl_gb` `nucl_wgs` from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/. Previously downloaded files can be directly provided with this argument.

When `--input-target file`, `--ncbi-file-info` uses `assembly_summary.txt` from https://ftp.ncbi.nlm.nih.gov/genomes/ to extract target information (options `refseq` `genbank` `refseq_historical` `genbank_historical`. Previously downloaded files can be directly provided with this argument.

If you are using outdated, removed or inactive assembly or sequence files and accessions from NCBI, make sure to include `dead_nucl` `dead_wgs` for `--ncbi-sequence-info` or `refseq_historical` `genbank_historical` for `--ncbi-file-info`. `eutils` option does not work with outdated accessions.

### ganon report

#### report type (--report-type)

Several reports are availble with `--report-type`: `reads`, `abundance`, `dist`, `corr`, `matches`:

`reads` reports **sequence abundances** which are the basic proportion of reads classified in the sample.

`abundance` will convert sequence abundance into **taxonomic abundances** by re-distributing read counts among leaf nodes and correcting by genome size. The re-distribution applies for reads classified with a LCA assignment and it is proportional to the number of unique matches of leaf nodes available in the ganon database (relative to the LCA node). Genome size is estimated based on [NCBI or GTDB auxiliary files](#genome-sizes---genome-size-files). Genome size correction is applied by rank based on default ranks only (superkingdom phylum class order family genus species assembly). Read counts in intermediate ranks will be corrected based on the closest parent default rank and re-assigned to its original rank.

`dist` is the same of `reads` with read count re-distribution

`corr` is the same of `reads` with correction by genome size

`matches` will report the total number of matches classified, either unique or shared. *This report will output the total number of matches instead the total number of reads reported in all other reports.*

## Parameters

```
usage: ganon [-h] [-v] {build,build-custom,update,classify,report,table} ...

- - - - - - - - - -
   _  _  _  _  _   
  (_|(_|| |(_)| |  
   _|   v. 1.3.0
- - - - - - - - - -

positional arguments:
  {build,build-custom,update,classify,report,table}
    build               Download and build ganon default databases
                        (refseq/genbank)
    build-custom        Build custom ganon databases
    update              Update ganon default databases
    classify            Classify reads against built databases
    report              Generate reports from classification results
    table               Generate table from reports

options:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
```

<details>
  <summary>ganon build</summary>

```
usage: ganon build [-h] [-g [...]] [-a [...]] [-b [...]] [-o] [-c] [-u] [-m [...]] [-z [...]] -d DB_PREFIX [-x] [-t]
                   [-p] [-f] [-k] [-w] [-s] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -g [ ...], --organism-group [ ...]
                        One or more organism groups to download [archaea,bacteria,fungi,human,invertebrate,metagenomes,o
                        ther,plant,protozoa,vertebrate_mammalian,vertebrate_other,viral]. Mutually exclusive --taxid
                        (default: None)
  -a [ ...], --taxid [ ...]
                        One or more taxonomic identifiers to download. e.g. 562 (-x ncbi) or 's__Escherichia coli' (-x
                        gtdb). Mutually exclusive --organism-group (default: None)
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (default: None)

download arguments:
  -b [ ...], --source [ ...]
                        Source to download [refseq,genbank] (default: ['refseq'])
  -o , --top            Download limited assemblies for each taxa. 0 for all. (default: 0)
  -c, --complete-genomes
                        Download only sub-set of complete genomes (default: False)
  -u , --genome-updater 
                        Additional genome_updater parameters (https://github.com/pirovc/genome_updater) (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)

important arguments:
  -x , --taxonomy       Set taxonomy to enable taxonomic classification, lca and reports [ncbi,gtdb,skip] (default:
                        ncbi)
  -t , --threads 

advanced arguments:
  -p , --max-fp         Max. false positive rate for bloom filters Mutually exclusive --filter-size. (default: 0.05)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. (default: 0)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [0-5]. 0 to detect optimal value.
                        (default: 4)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon build-custom</summary>

```
usage: ganon build-custom [-h] [-i [...]] [-e] [-n] [-a] [-l] [-m [...]] [-z [...]] [-r [...]] [-q [...]] -d DB_PREFIX
                          [-x] [-t] [-p] [-f] [-k] [-w] [-s] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). Mutually exclusive --input-file. (default: None)
  -e , --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: fna.gz)
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (default: None)

custom arguments:
  -n , --input-file     Manually set information for input files: file <tab> [target <tab> node <tab> specialization
                        <tab> specialization name]. target is the sequence identifier if --input-target sequence (file
                        can be repeated for multiple sequences). if --input-target file and target is not set, filename
                        is used. node is the taxonomic identifier. Mutually exclusive --input (default: None)
  -a , --input-target   Target to use [file, sequence]. By default: 'file' if multiple input files are provided or
                        --input-file is set, 'sequence' if a single file is provided. Using 'file' is recommended and
                        will speed-up the building process (default: None)
  -l , --level          Use a specialized target to build the database. By default, --level is the --input-target.
                        Options: any available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy).
                        Further specialization options [assembly,custom]. assembly will retrieve and use the assembly
                        accession and name. custom requires and uses the specialization field in the --input-file.
                        (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)

ncbi arguments:
  -r [ ...], --ncbi-sequence-info [ ...]
                        Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information.
                        [eutils,nucl_gb,nucl_wgs,nucl_est,nucl_gss,pdb,prot,dead_nucl,dead_wgs,dead_prot or one or more
                        accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/]. By
                        default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise. (default:
                        [])
  -q [ ...], --ncbi-file-info [ ...]
                        Downloads assembly_summary files to extract target information.
                        [refseq,genbank,refseq_historical,genbank_historical or one or more assembly_summary files from
                        https://ftp.ncbi.nlm.nih.gov/genomes/] (default: ['refseq', 'genbank'])

important arguments:
  -x , --taxonomy       Set taxonomy to enable taxonomic classification, lca and reports [ncbi,gtdb,skip] (default:
                        ncbi)
  -t , --threads 

advanced arguments:
  -p , --max-fp         Max. false positive rate for bloom filters Mutually exclusive --filter-size. (default: 0.05)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. (default: 0)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [0-5]. 0 to detect optimal value.
                        (default: 4)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon update</summary>

```
usage: ganon update [-h] -d DB_PREFIX [-o] [-t] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Existing database input prefix (default: None)

important arguments:
  -o , --output-db-prefix 
                        Output database prefix. By default will be the same as --db-prefix and overwrite files (default:
                        None)
  -t , --threads 

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon classify</summary>

```
usage: ganon classify [-h] -d [DB_PREFIX ...] [-s [reads.fq[.gz] ...]] [-p [reads.1.fq[.gz] reads.2.fq[.gz] ...]]
                      [-c [...]] [-e [...]] [-o] [--output-lca] [--output-all] [--output-unclassified] [--output-single]
                      [-t] [-l [...]] [-r [...]] [--verbose] [--quiet]

options:
  -h, --help            show this help message and exit

required arguments:
  -d [DB_PREFIX ...], --db-prefix [DB_PREFIX ...]
                        Database input prefix[es] (default: None)
  -s [reads.fq[.gz] ...], --single-reads [reads.fq[.gz] ...]
                        Multi-fastq[.gz] file[s] to classify (default: None)
  -p [reads.1.fq[.gz] reads.2.fq[.gz] ...], --paired-reads [reads.1.fq[.gz] reads.2.fq[.gz] ...]
                        Multi-fastq[.gz] pairs of file[s] to classify (default: None)

cutoff/filter arguments:
  -c [ ...], --rel-cutoff [ ...]
                        Min. percentage of a read (set of minimizers) shared with the a reference necessary to consider
                        a match. Generally used to cutoff low similarity matches. Single value or one per database (e.g.
                        0.7 1 0.25). 0 for no cutoff (default: [0.75])
  -e [ ...], --rel-filter [ ...]
                        Additional relative percentage of minimizers (relative to the best match) to keep a match.
                        Generally used to select best matches above cutoff. Single value or one per hierarchy (e.g. 0.1
                        0). 1 for no filter (default: [0.0])

output arguments:
  -o , --output-prefix 
                        Output prefix for output (.rep) and report (.tre). Empty to output to STDOUT (only .rep)
                        (default: None)
  --output-lca          Output an additional file with one lca match for each read (.lca) (default: False)
  --output-all          Output an additional file with all matches. File can be very large (.all) (default: False)
  --output-unclassified
                        Output an additional file with unclassified read headers (.unc) (default: False)
  --output-single       When using multiple hierarchical levels, output everything in one file instead of one per
                        hierarchy (default: False)

other arguments:
  -t , --threads        Number of sub-processes/threads to use (default: 1)
  -l [ ...], --hierarchy-labels [ ...]
                        Hierarchy definition of --db-prefix files to be classified. Can also be a string, but input will
                        be sorted to define order (e.g. 1 1 2 3). The default value reported without hierarchy is 'H1'
                        (default: None)
  -r [ ...], --ranks [ ...]
                        Ranks to report taxonomic abundances (.tre). empty will report default ranks
                        [superkingdom,phylum,class,order,family,genus,species,assembly]. This file can be re-generated
                        with the 'ganon report' command for other types of abundances (reads, matches) with further
                        filtration and output options (default: [])
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
```

</details>

<details>
  <summary>ganon report</summary>

```
usage: ganon report [-h] -i [...] [-e INPUT_EXTENSION] -o OUTPUT_PREFIX [-d [...]] [-x] [-m [...]] [-z [...]] [-f] [-t]
                    [-r [...]] [-s] [-a] [-y] [-p [...]] [-k [...]] [-c] [--verbose] [--quiet] [--min-count]
                    [--max-count] [--names [...]] [--names-with [...]] [--taxids [...]]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). '.rep' file(s) from ganon classify. (default: None)
  -e INPUT_EXTENSION, --input-extension INPUT_EXTENSION
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: rep)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix for report file 'output_prefix.tre'. In case of multiple files, the base input
                        filename will be appended at the end of the output file 'output_prefix + FILENAME.tre' (default:
                        None)

db/tax arguments:
  -d [ ...], --db-prefix [ ...]
                        Database prefix(es) used for classification. Only '.tax' file(s) are required. If not provided,
                        new taxonomy will be downloaded. Mutually exclusive with --taxonomy. (default: [])
  -x , --taxonomy       Taxonomy database to use [ncbi,gtdb,skip]. Mutually exclusive with --db-prefix. (default: ncbi)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)

output arguments:
  -f , --output-format 
                        Output format [text,tsv,csv,bioboxes]. text outputs a tabulated formatted text file for better
                        visualization. bioboxes is the the CAMI challenge profiling format (only percentage/abundances
                        are reported). (default: tsv)
  -t , --report-type    Type of report [abundance,reads,matches,dist,corr]. 'abundance' -> tax. abundance (re-distribute
                        read counts and correct by genome size), 'reads' -> sequence abundance, 'matches' -> report all
                        unique and shared matches, 'dist' -> like reads with re-distribution of shared read counts only,
                        'corr' -> like abundance without re-distribution of shared read counts (default: abundance)
  -r [ ...], --ranks [ ...]
                        Ranks to report ['', 'all', custom list]. 'all' for all possible ranks. empty for default ranks
                        [superkingdom,phylum,class,order,family,genus,species,assembly]. (default: [])
  -s , --sort           Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage
                        (with --ranks all) (default: )
  -a, --no-orphan       Omit orphan nodes from the final report. Otherwise, orphan nodes (= nodes not found in the
                        db/tax) are reported as 'na' with root as direct parent. (default: False)
  -y, --split-hierarchy
                        Split output reports by hierarchy (from ganon classify --hierarchy-labels). If activated, the
                        output files will be named as '{output_prefix}.{hierarchy}.tre' (default: False)
  -p [ ...], --skip-hierarchy [ ...]
                        One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels) (default:
                        [])
  -k [ ...], --keep-hierarchy [ ...]
                        One or more hierarchies to keep in the report (from ganon classify --hierarchy-labels) (default:
                        [])
  -c , --top-percentile 
                        Top percentile filter, based on percentage/relative abundance. Applied only at default ranks
                        [superkingdom,phylum,class,order,family,genus,species,assembly] (default: 0)

optional arguments:
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)

filter arguments:
  --min-count           Minimum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --max-count           Maximum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --names [ ...]        Show only entries matching exact names of the provided list (default: [])
  --names-with [ ...]   Show entries containing full or partial names of the provided list (default: [])
  --taxids [ ...]       One or more taxids to report (including children taxa) (default: [])
```

</details>

<details>
  <summary>ganon table</summary>

```
usage: ganon table [-h] -i [...] [-e] -o OUTPUT_FILE [-l] [-f] [-t] [-a] [-m] [-r] [-n] [--header]
                   [--unclassified-label] [--filtered-label] [--skip-zeros] [--transpose] [--verbose] [--quiet]
                   [--min-count] [--max-count] [--names [...]] [--names-with [...]] [--taxids [...]]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). '.tre' file(s) from ganon report. (default: None)
  -e , --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: tre)
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output filename for the table (default: None)

output arguments:
  -l , --output-value   Output value on the table [percentage, counts]. percentage values are reported between [0-1]
                        (default: counts)
  -f , --output-format 
                        Output format [tsv, csv] (default: tsv)
  -t , --top-sample     Top hits of each sample individually (default: 0)
  -a , --top-all        Top hits of all samples (ranked by percentage) (default: 0)
  -m , --min-frequency 
                        Minimum number/percentage of files containing an taxa to keep the taxa [values between 0-1 for
                        percentage, >1 specific number] (default: 0)
  -r , --rank           Define specific rank to report. Empty will report all ranks. (default: None)
  -n, --no-root         Do not report root node entry and lineage. Direct and shared matches to root will be accounted
                        as unclassified (default: False)
  --header              Header information [name, taxid, lineage] (default: name)
  --unclassified-label 
                        Add column with unclassified count/percentage with the chosen label. May be the same as
                        --filtered-label (e.g. unassigned) (default: None)
  --filtered-label      Add column with filtered count/percentage with the chosen label. May be the same as
                        --unclassified-label (e.g. unassigned) (default: None)
  --skip-zeros          Do not print lines with only zero count/percentage (default: False)
  --transpose           Transpose output table (taxa as cols and files as rows) (default: False)

optional arguments:
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)

filter arguments:
  --min-count           Minimum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --max-count           Maximum number/percentage of counts to keep an taxa [values between 0-1 for percentage, >1
                        specific number] (default: 0)
  --names [ ...]        Show only entries matching exact names of the provided list (default: [])
  --names-with [ ...]   Show entries containing full or partial names of the provided list (default: [])
  --taxids [ ...]       One or more taxids to report (including children taxa) (default: [])
```

</details>
