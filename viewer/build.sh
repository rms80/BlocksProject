#!/bin/bash
set -e
cd "$(dirname "$0")"

# Configure if needed
if [ ! -f build/Makefile ]; then
    cmake -S . -B build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
fi

# Build
cmake --build build --target polyscope_viewer -j$(sysctl -n hw.ncpu)
