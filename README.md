# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon)

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
