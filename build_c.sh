#!/usr/bin/env bash
set -e

rm -rf util_c/build
mkdir -p util_c/build
cd util_c/build
cmake ..
make

cd ../..
rm -rf util_c/bin util_c/lib
mkdir util_c/bin util_c/lib

mv util_c/build/bin/* util_c/bin/
mv util_c/build/lib/* util_c/lib/

rm -rf util_c/build

echo "✅ Binaries in ./bin/, shared libraries in ./lib/"