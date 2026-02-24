#!/bin/bash
set -e
cd "$(dirname "$0")"

# Configure if needed
if [ ! -f build/Makefile ]; then
    cmake -S . -B build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
fi

# Build cell_match_test
cmake --build build --target cell_match_test -j$(sysctl -n hw.ncpu)

# Run from build dir so ../output/ resolves to BlocksProject/output/
mkdir -p output
echo "--- Running cell_match_test ---"
cd build
./cell_match_test
