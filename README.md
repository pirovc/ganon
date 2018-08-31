## Dependencies

- gcc7
- python3.5
- cmake3

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
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-std=c++1z -Wuninitialized -W -Wall -Wstrict-aliasing -pedantic -Wno-long-long -Wno-variadic-macros -Wunused -msse4.2 -O3 -DNDEBUG -static -pthread -march=native" -DSEQAN_INCLUDE_PATH="${SEQAN_PATH}/include" -DSDSL_INCLUDE_DIRS="${SDSL_PATH}/include" -DCMAKE_PREFIX_PATH="${SEQAN_PATH}/util/cmake" ../

	make

## Usage

	./ganon -h
	./ganon build -h
	./ganon classify -h
	./ganon update -h

## Installing GCC7 in a separate environment with conda

	conda create -n gcc7 -c quantstack gcc-7 libgcc-7
	source activate gcc7
