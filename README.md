# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon)

ganon is a k-mer based read classification tool which uses Interleaved Bloom Filters in conjunction with a taxonomic clustering and a k-mer counting-filtering scheme. 

> **ganon: continuously up-to-date with database growth for precise short read classification in metagenomics**
> Vitor C. Piro, Temesgen H. Dadi, Enrico Seiler, Knut Reinert, Bernhard Y. Renard
> bioRxiv 406017; doi: [10.1101/406017](https://doi.org/10.1101/406017)

## Installation

```shh
conda install -c bioconda ganon
```

To install ganon directly from the source, please check [manual installation](#manual-installation))

## Running ganon with sample data

### build

```shh
ganon build --db-prefix sample_bacteria --input-files tests/ganon-build/data/sequences/bacteria*.fasta.gz
```

### classify

```shh
ganon classify --db-prefix sample_bacteria --reads tests/ganon-classify/data/reads/bacteria.simulated.1.fq -o sample_results
```

### update

```
ganon update --db-prefix sample_bacteria --output-db-prefix sample_bateria_virus --input-files tests/ganon-build/data/sequences/virus*.fasta.gz
```

## Output files

The main output file is the `*.tre` which will sumarize the results.

Example of the sample data classification `sample_results.tre` previously generated:

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

running with `--ranks all`, the output will show all ranks used for classification, sorted by lineage:

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

### classify

#### .out

all matches (one per hierarchy)

	readid <tab> assignment <tab> k-mer count

* a negative k-mer count work a flag to show that such read did not have the minimum amount of matches to pass the `--max-error-unique` threshold

#### .lca

only one match / read (one per hierarchy)
	
	readid <tab> lca assignment <tab> max k-mer count

LCA script obtained from https://www.ics.uci.edu/~eppstein/

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
taxonomic information required for LCA and reports

## Manual Installation

### Dependencies

#### build

- gcc >=7 (check [gcc7 with conda](#installing-gcc7-in-a-separate-environment-with-conda))
- cmake >=3
- Catch2 >=2.7.0
- cxxopts >=2.1.2
- sdsl-lite 3.0 [d6ed14](https://github.com/xxsds/sdsl-lite/commit/d6ed14d5d731ed4a4ec12627c1ed7154b396af48)
- seqan 2.4.0 [c308e9](https://github.com/eseiler/seqan/commit/c308e99f10d942382d4c7ed6fc91be1a889e644c)

#### run

- python >=3.5
- taxsbp >=0.1.1
- binpacking ==1.4.1
- pandas

### Cloning ganon and taxsbp

```shh
git clone --recurse-submodules https://github.com/pirovc/ganon.git
git clone https://github.com/pirovc/taxsbp.git
```

To compile a specific ganon release, check [Using a specific release](#using-a-specific-release))

### Installing taxsbp + binpacking

```shh
cd taxsbp
python setup.py install
taxsbp -h
```

### Building (ganon-build and ganon-classify)
	
```shh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Installing GCC7 in a separate environment with conda

```shh
conda create -n gcc7 -c quantstack gcc-7 libgcc-7 cmake>=3.8.2
source activate gcc7
```

If you are getting the following error `ganon-classify: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'CXXABI_1.3.11' not found` you have to set the `LD_LIBRARY_PATH`

```shh
export LD_LIBRARY_PATH=/home/user/miniconda3/envs/gcc7/lib/
```
## Using a specific release

```shh
# download and unpack release
wget https://github.com/pirovc/ganon/archive/0.1.0.tar.gz
tar xf 0.1.0.tar.gz
cd ganon-0.1.0

# get submodules
rm -r libs/*
git init
git config -f .gitmodules --get-regexp '^submodule\..*\.path$' | 
  while read path_key path; do
    url=$(git config -f .gitmodules --get "$(echo $path_key | sed 's/\.path/.url/')")
    branch=$(git config -f .gitmodules --get "$(echo $path_key | sed 's/\.path/.branch/')")
    if [[ ! -z $branch ]]; then branch="-b ${branch}"; fi
    git submodule add $branch $url $path
  done
```

Make sure you commit the sdsl-lite and seqan libs in the required commits