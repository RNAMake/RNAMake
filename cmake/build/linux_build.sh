#!/bin/sh
set -e
# delete the old stuff
echo "Deleting things from previous builds..."

rm -f CMakeCache.txt
rm -rf -f CMakeFiles/

echo "Rebuilding CMakeLists.txt file..."

python ../make_project.py -target "linux" # this should have a mac flag

echo "Calling cmake..."

cmake .

echo "Assembling in make..."

make

echo "Moving binaries to $RNAMAKE/linux-binaries/"

mv *unittest ../../linux-binaries/

mv ddg_tecto design_rna_scaffold eternabot get_best_solutions sequence_optimization_benchmarks sequence_optimizer_app thermo_simulation ../../linux-binaries/


echo "Build complete!"
