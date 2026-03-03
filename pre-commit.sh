#!/bin/bash
ruff format
ruff check --fix

# C++
find src/ -iname '*.hpp' -o -iname '*.cpp' | xargs clang-format -i