#!/bin/bash
 # note you will have to escape special characters if they are in your password.
git clone --branch [BRANCH] https://[USERNAME]:[PASSWORD]@github.com/RNAMake/RNAMake.git && \
    git checkout jdy/update_build_motif_graph \
    export RNAMAKE=/RNAMake/ && \
    cd RNAMake/cmake/build/ && \
    python ../make_project.py -target mac && \
    cmake -G Ninja
    time ninja
    cd /
    # rm -rf -f *.a CMakeCache.txt CMakeFiles/ *.txt *.cmake ../../src/ ../../.git* && \
    #cd / && \
    tar -cvzf RNAMake-binaries.tar.gz RNAMake/
