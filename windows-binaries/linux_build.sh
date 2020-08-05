#!/bin/sh

set -e
mkdir -p /result
# delete the old stuff
echo "Deleting things from previous builds..."

git    clone  --branch client-binaries https://chrisjurich:Ilikesoccer5!@github.com/jyesselm/RNAMake.git

#rm -f CMakeCache.txt
#rm -rf -f CMakeFiles/
#

pushd RNAMake//cmake/build
#echo "Rebuilding CMakeLists.txt file..."
#
python ../make_project.py -target "linux" # this should have a mac flag
#
#echo "Calling cmake..."
#
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .
#
#echo "Assembling in make..."
#
time make
#
#echo "Moving binaries to $RNAMAKE/linux-binaries/"
#
#mv *unittest ../../linux-binaries/
#
#mv ddg_tecto design_rna_scaffold eternabot get_best_solutions sequence_optimization_benchmarks sequence_optimizer_app thermo_simulation ../../linux-binaries/
#
#
#echo "Build complete!"
#!/bin/bash

rm -rf CMakeFiles/
pushd /RNAMake//
rm -rf src/
pushd /
tar -czvf rnamake-binaries-linux.tar.gz RNAMake/cmake/build/
