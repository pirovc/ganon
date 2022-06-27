# ganon [![Build Status](https://travis-ci.com/pirovc/ganon.svg?branch=master)](https://travis-ci.com/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

ganon classifies short reads against large sets of refence sequences efficiently, with custom download and build, taxonomic classification (ncbi and gtdb) and many other [features](#Features). [Publication 10.1101/406017](https://dx.doi.org/10.1093/bioinformatics/btaa458)

## Quick install/usage guide

#### Install with conda

```sh
conda install -c bioconda -c conda-forge ganon
```

#### Download and build
```sh
# Archaeal complete genome sequences from NCBI RefSeq
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 12
```

#### Classify
```sh
ganon classify --db-prefix arc_cg_rs --output-prefix classify_results --single-reads my_reads.fq.gz --threads 12
```

#### Generate filtered reports
```sh
ganon report --db-prefix arc_cg_rs --input classify_results.rep --output-prefix filtered_report --min-count 0.01
```

## Details

ganon is designed to index large sets of genomic reference sequences and to classify short reads against them efficiently. The tool uses Interleaved Bloom Filters as indices based on k-mer and minimizer sequences. It was mainly developed, but not limited, to the metagenomic classification problem: assign short fragments to their closest reference among thousands of references.

### Features

- NCBI and GTDB native support for taxonomic classification
- integrated download of commonly used reference sequences from RefSeq/Genbank (`ganon build`)
- [update indices](#updating-the-index) incrementally (`ganon update`)
- customizable build for non-standard sequence files (`ganon build-custom`)
- build and classify at different taxonomic levels, strain/assembly or custom specialization
- perform [hierarchical classification](#multiple-and-hierarchical-classification)
- report the lowest common ancestor (LCA), multiple and unique matches for every read
- generate reports and tables for multi-sample studies with filtering options

## Installation guide

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

```bash
conda install -c bioconda -c conda-forge ganon
```

* There are possible performance benefits compiling ganon from source in the target machine rather than using the conda version. To do so, please follow the instructions below:

<details>
  <summary>Instructions</summary>

### build dependencies

System packages:
- gcc >=8
- cmake >=3.10
- zlib

### run dependencies

System packages:
- python >=3.5
- pandas >=0.22.0
- multitax >=1.1.1

```bash
python3 -m pip install "pandas>=0.22.0"
python3 -m pip install "multitax==1.1.1"
```

### Downloading and building ganon

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

### Run tests

```bash
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/
cd build_cpp/
ctest -VV .
```

</details>

## Usage guide

save states

### Examples

### Building custom indices

### Updating indices

## Output files

### build/update

Every run on `ganon build` or `ganon update` will generate the following database files:

 - {prefix}**.ibf**: main interleaved bloom filter index file
 - {prefix}**.tax**: taxonomic tree *(fields: target/node, parent, rank, name)* (if --taxonomy is used)

Obs:
-  Database files generated with version 1.2.0 or higher are not compatible with older versions.

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

</details>

## Multiple and Hierarchical classification

## Choosing/explaining parameters

## Parameters