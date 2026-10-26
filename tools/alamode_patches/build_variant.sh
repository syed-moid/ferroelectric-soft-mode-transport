#!/bin/sh
# Build anphon in one or more ALAMODE trees with the conda-forge toolchain of env.explicit.txt.
# Usage: TOOLCHAIN=<dir holding env/ and the trees> sh build_variant.sh <tree> [<tree> ...]
# The recipe is the one used for every build listed in BUILD.md (only the tree name differs).
set -e
T=${TOOLCHAIN:?set TOOLCHAIN to the directory that holds env/ and the ALAMODE trees}
export PATH=$T/env/bin:$PATH
cd "$T"
for tree in "$@"; do
cmake -S $tree -B $tree/_build -DCMAKE_BUILD_TYPE=Release  -DCMAKE_PREFIX_PATH=$T/env -DCMAKE_C_COMPILER=$T/env/bin/clang -DCMAKE_CXX_COMPILER=$T/env/bin/clang++  -DCMAKE_CXX_FLAGS="-I$T/env/include/eigen3 -I$T/env/include -fopenmp" -DCMAKE_C_FLAGS="-I$T/env/include -fopenmp"  -DCMAKE_EXE_LINKER_FLAGS="-L$T/env/lib -lomp -Wl,-rpath,$T/env/lib" -DSPGLIB_ROOT=$T/env -DFFTW3_ROOT=$T/env
cmake --build $tree/_build -j12 --target anphon
done
