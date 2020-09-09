#!/bin/bash
 # note you will have to escape special characters if they are in your password.
git clone --branch [BRANCH] https://[USERNAME]:[PASSWORD]@github.com/RNAMake/RNAMake.git && \
    export RNAMAKE=/RNAMake/ && \
    export X3DNA=$RNAMAKE/resources/x3dna/linux/ && \
    cd RNAMake/cmake/build/ && \
    python ../make_project.py -target mac && \
    cmake -G Ninja
    time ninja
    cd /
    # rm -rf -f *.a CMakeCache.txt CMakeFiles/ *.txt *.cmake ../../src/ ../../.git* && \
    #cd / && \
    tar -cvzf RNAMake-binaries.tar.gz RNAMake/
