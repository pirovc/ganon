#!/bin/bash
echo "---- ruff format/check ----"
ruff format
ruff check --fix

# C++
echo "---- clang-format ----"
find {src/,tests/} -iname '*.hpp' -o -iname '*.cpp' | xargs clang-format --verbose -i

# Update docs
echo "---- update docs ----"
docs/generate_params.sh > docs/params.md
v=$(ganon --version | cut -d" " -f 3)
sed -i 's/'$v'/'${v:0:5}'/' docs/params.md

# Tests
echo "---- tests c++ ----"
ctest --test-dir build/
echo "---- tests python ----"
python3 -m unittest discover -s tests/ganon/integration/