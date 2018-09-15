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

	# taxsbp
	pip install binpacking==1.3
	git clone https://github.com/pirovc/taxsbp/

	# seqan
	git clone https://github.com/xxsds/sdsl-lite/
	SDSL_PATH=$(readlink -e sdsl-lite)
	git clone -b unified https://github.com/eseiler/seqan.git
	SEQAN_PATH=$(readlink -e seqan)
	
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make

## Usage

	./ganon -h
	./ganon build -h
	./ganon classify -h
	./ganon update -h

## Installing GCC7 in a separate environment with conda

	conda create -n gcc7 -c quantstack gcc-7 libgcc-7
	source activate gcc7
