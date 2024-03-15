#!/usr/bin/bash

rm -r build
mkdir build
cd build
g++ ../MonteCarloEfficiency.cxx -o MonteCarloEfficiency -std=c++11 -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I ../../star-upc-new/include/ -I ../include -L ../../star-upc-new/build/ -lstar-upc
export LD_LIBRARY_PATH=../star-upc-new/build/:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=../star-upc-new/build/:${DYLD_LIBRARY_PATH}
