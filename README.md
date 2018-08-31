Dependencies

gcc7
python3.5
cmake3

Installation

	# taxsbp
	pip install binpacking==1.3
	git clone https://github.com/pirovc/taxsbp/

	# seqan
	git clone https://github.com/xxsds/sdsl-lite/
	git clone --single-branch -b unified https://github.com/eseiler/seqan.git
	
	mkdir build
	cd build

	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-std=c++1z -Wuninitialized -W -Wall -Wstrict-aliasing -pedantic -Wno-long-long -Wno-variadic-macros -Wunused -msse4.2 -O3 -DNDEBUG -static -pthread -march=native" -DSEQAN_INCLUDE_PATH="../seqan/include" -DSDSL_INCLUDE_DIRS="../sdsl-lite/include" -DCMAKE_PREFIX_PATH="../seqan/util/cmake" ../

	make

Usage

		