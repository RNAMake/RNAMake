# using Ubuntu for the base image
FROM ubuntu:latest
FROM python:3.7.7

MAINTAINER Chris Jurich <cjurich2@huskers.unl.edu>

RUN apt-get update && 				\
	apt-get install -y  mingw-w64 		\
	cmake 					\
	file 					\
	libc6-dev


COPY . /
RUN git  clone    --branch client-binaries https://chrisjurich:Ilikesoccer5!@github.com/jyesselm/RNAMake.git
WORKDIR /RNAMake/cmake/build/ 
RUN python ../make_project.py -target "windows" &&                      \
    cmake . -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc  &&          \
    make 
RUN rm -rf CMakeFiles/                          &&                      \
    cd /RNAMake//                            &&                      \
    rm -rf src/                                 &&                      \
    cd /                                     &&                      \
    tar -czvf windows-binaries.tar.gz RNAMake/cmake/build/

ENTRYPOINT cat windows-binaries.tar.gz

#RUN python ../make_project.py
#RUN cmake .
#RUN pwd 
#RUN make

#set -e
#mkdir -p /result
## delete the old stuff
#echo "Deleting things from previous builds..."
#
#
##rm -f CMakeCache.txt
##rm -rf -f CMakeFiles/
##
#
#pushd RNAMake//cmake/build
##echo "Rebuilding CMakeLists.txt file..."
##
#python ../make_project.py -target "linux" # this should have a mac flag
##
##echo "Calling cmake..."
##
#cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .
##
##echo "Assembling in make..."
##
#time make
##
##echo "Moving binaries to $RNAMAKE/linux-binaries/"
##
##mv *unittest ../../linux-binaries/
##
##mv ddg_tecto design_rna_scaffold eternabot get_best_solutions sequence_optimization_benchmarks sequence_optimizer_app thermo_simulation ../../linux-binaries/
##
##
##echo "Build complete!"
##!/bin/bash
#
#rm -rf CMakeFiles/
#pushd /RNAMake//
#rm -rf src/
#pushd /
#tar -czvf rnamake-binaries-linux.tar.gz RNAMake/cmake/build/
