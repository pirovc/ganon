# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon)


## Dependencies

- gcc7 (check [gcc7 with conda](#installing-gcc7-in-a-separate-environment-with-conda))
- python3.5
- cmake3

## Cloning

Make sure to clone the repository with its submodules. One way to do this is as follows:

```shh
git clone --recurse-submodules https://github.com/pirovc/ganon.git
```

## Installation

Installing *binpacking* and *taxsbp*:

```shh
pip install binpacking==1.3
git clone https://github.com/pirovc/taxsbp/
```

## Building
	
```shh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Runnin ganon with sample data

### build

```shh
./ganon build --db-prefix sample_bacteria --input-files build/data/sequences/bacteria*.fasta.gz     # --ganon-path build/ --taxsbp-path taxsbp/
```

* it may be necessary to set-up the path for the ganon binaries `--ganon-path` and taxsbp `--taxsbp-path`

* ganon will automatically retrieve taxonomic information from NCBI web services and ftp (it may take a while if you have too many sequences). To speed-up such process it is possible to manually provide the `--len-taxid-file` which is a tab separated file containing `sequence identifier <tab> sequence length <tab> taxonomic id [<tab> assembly/strain id]`

### classify

```shh
./ganon classify --db-prefix sample_bacteria --reads build/data/reads/bacteria.simulated.1.fq -o sample_results      # --ganon-path build/
```

`sample_results.tre` will sumarize the results (with fixed ranks):

```
unclassified  -        -                                              -                               1   1.00000
root          1        1                                              root                            99  99.00000
superkingdom  2        1|2                                            Bacteria                        99  99.00000
phylum        1239     1|2|1239                                       Firmicutes                      57  57.00000
phylum        1224     1|2|1224                                       Proteobacteria                  42  42.00000
class         91061    1|2|1239|91061                                 Bacilli                         57  57.00000
class         28211    1|2|1224|28211                                 Alphaproteobacteria             29  29.00000
class         1236     1|2|1224|1236                                  Gammaproteobacteria             13  13.00000
order         1385     1|2|1239|91061|1385                            Bacillales                      57  57.00000
order         204458   1|2|1224|28211|204458                          Caulobacterales                 29  29.00000
order         72274    1|2|1224|1236|72274                            Pseudomonadales                 13  13.00000
family        186822   1|2|1239|91061|1385|186822                     Paenibacillaceae                57  57.00000
family        76892    1|2|1224|28211|204458|76892                    Caulobacteraceae                29  29.00000
family        468      1|2|1224|1236|72274|468                        Moraxellaceae                   13  13.00000
genus         44249    1|2|1239|91061|1385|186822|44249               Paenibacillus                   57  57.00000
genus         75       1|2|1224|28211|204458|76892|75                 Caulobacter                     29  29.00000
genus         469      1|2|1224|1236|72274|468|469                    Acinetobacter                   13  13.00000
species       1406     1|2|1239|91061|1385|186822|44249|1406          Paenibacillus polymyxa          57  57.00000
species       366602   1|2|1224|28211|204458|76892|75|366602          Caulobacter sp. K31             29  29.00000
species       470      1|2|1224|1236|72274|468|469|470                Acinetobacter baumannii         13  13.00000
species+      1052684  1|2|1239|91061|1385|186822|44249|1406|1052684  Paenibacillus polymyxa M1       57  57.00000
species+      696749   1|2|1224|1236|72274|468|469|470|696749         Acinetobacter baumannii 1656-2  13  13.00000
```

if you run the classification with `--ranks all`, the output will show all ranks used for classification, sorted by lineage:

```
unclassified   -        -                                                             -                                              1   1.00000
no rank        1        1                                                             root                                           99  99.00000
no rank        131567   1|131567                                                      cellular organisms                             99  99.00000
superkingdom   2        1|131567|2                                                    Bacteria                                       99  99.00000
phylum         1224     1|131567|2|1224                                               Proteobacteria                                 42  42.00000
class          1236     1|131567|2|1224|1236                                          Gammaproteobacteria                            13  13.00000
order          72274    1|131567|2|1224|1236|72274                                    Pseudomonadales                                13  13.00000
family         468      1|131567|2|1224|1236|72274|468                                Moraxellaceae                                  13  13.00000
genus          469      1|131567|2|1224|1236|72274|468|469                            Acinetobacter                                  13  13.00000
species group  909768   1|131567|2|1224|1236|72274|468|469|909768                     Acinetobacter calcoaceticus/baumannii complex  13  13.00000
species        470      1|131567|2|1224|1236|72274|468|469|909768|470                 Acinetobacter baumannii                        13  13.00000
no rank        696749   1|131567|2|1224|1236|72274|468|469|909768|470|696749          Acinetobacter baumannii 1656-2                 13  13.00000
class          28211    1|131567|2|1224|28211                                         Alphaproteobacteria                            29  29.00000
order          204458   1|131567|2|1224|28211|204458                                  Caulobacterales                                29  29.00000
family         76892    1|131567|2|1224|28211|204458|76892                            Caulobacteraceae                               29  29.00000
genus          75       1|131567|2|1224|28211|204458|76892|75                         Caulobacter                                    29  29.00000
species        366602   1|131567|2|1224|28211|204458|76892|75|366602                  Caulobacter sp. K31                            29  29.00000
no rank        1783272  1|131567|2|1783272                                            Terrabacteria group                            57  57.00000
phylum         1239     1|131567|2|1783272|1239                                       Firmicutes                                     57  57.00000
class          91061    1|131567|2|1783272|1239|91061                                 Bacilli                                        57  57.00000
order          1385     1|131567|2|1783272|1239|91061|1385                            Bacillales                                     57  57.00000
family         186822   1|131567|2|1783272|1239|91061|1385|186822                     Paenibacillaceae                               57  57.00000
genus          44249    1|131567|2|1783272|1239|91061|1385|186822|44249               Paenibacillus                                  57  57.00000
species        1406     1|131567|2|1783272|1239|91061|1385|186822|44249|1406          Paenibacillus polymyxa                         57  57.00000
no rank        1052684  1|131567|2|1783272|1239|91061|1385|186822|44249|1406|1052684  Paenibacillus polymyxa M1                      57  57.00000
```

### update

```
./ganon update --db-prefix sample_bacteria --output-db-prefix sample_bateria_virus --input-files build/data/sequences/virus*.fasta.gz        # --ganon-path build/ --taxsbp-path taxsbp/
```

## Output files

### classify

#### .out

all matches (one per hierarchy)

	readid <tab> assignment <tab> k-mer count

* a negative k-mer count work a flag to show that such read did not have the minimum amount of matches to pass the `--max-error-unique` threshold

#### .lca

only one match / read (one per hierarchy)
	
	readid <tab> lca assignment <tab> max k-mer count

#### .rep

detailed report for lca matches (one per hierarchy)
	
	1) lca assignment <tab>
	2) reads assigned (lca) <tab>
	3) % reads assigned <tab>
	4) # reads assigned (total) <tab>
	5) # reads uniquely assigned <tab>
	6) sum of k-mers matched <tab>
	7) taxonomic rank <tab>
	8) name

#### .tre

tree-like output with cummulative counts and lineage (one per run)
	
	1) rank <tab>
	2) lca assignment <tab>
	3) lineage <tab>
	4) name <tab>
	5) cummulative # reads assigned <tab>
	6) cummulative % reads assigned

### build/update

#### .filter
main bloom filter
#### .map
maps assignment to binids
#### .bins
bins generated by TaxSBP
#### .nodes
taxonomic information required for LCA and LCA

## Installing GCC7 in a separate environment with conda

```shh
conda create -n gcc7 -c quantstack gcc-7 libgcc-7
source activate gcc7
```

If you are getting the following error `ganon-classify: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'CXXABI_1.3.11' not found` you have to set the `LD_LIBRARY_PATH`


```shh
export LD_LIBRARY_PATH=/home/user/miniconda3/envs/gcc7/lib/
```
