# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon)


## Dependencies

- gcc7
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
git clone -b v0.1 https://github.com/pirovc/taxsbp/
```

## Building
	
```shh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Output

ganon-classify:

.out (all matches)

	readid <tab> assignment <tab> k-mer count

.lca (one match / read)
	
	readid <tab> lca assignment <tab> max k-mer count

.rep (detailed report for lca matches)
	
	lca assignment <tab> # reads assigned (lca) <tab> % reads assigned <tab> # reads assigned (total) <tab> # reads uniquely assigned <tab> sum of k-mers matched <tab> taxonomic rank <tab> name

.tre (cummulative counts on the taxonomy tree)
	
	rank <tab> lca assignment <tab> lineage <tab> name <tab> cummulative # reads assigned <tab> cummulative % reads assigned

## Usage

```shh
./ganon -h
./ganon build -h
./ganon classify -h
./ganon update -h
```

## Installing GCC7 in a separate environment with conda

```shh
conda create -n gcc7 -c quantstack gcc-7 libgcc-7
source activate gcc7
```
