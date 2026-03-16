#!/bin/bash
ruff format
ruff check --fix

# C++
find src/ -iname '*.hpp' -o -iname '*.cpp' | xargs clang-format -i
find tests/ -iname '*.hpp' -o -iname '*.cpp' | xargs clang-format -i

# Update docs
docs/generate_params.sh > docs/params.md
v=$(ganon --version | cut -d" " -f 3)
sed -i 's/'$v'/'${v:0:5}'/' docs/params.md