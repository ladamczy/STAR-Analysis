#!/usr/bin/bash

rm -r build
mkdir build
cd build
g++ ../MonteCarloEfficiency.cxx -o MonteCarloEfficiency -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I ../../star-upc/include/ -I ../include -L ../../star-upc/build/ -lstar-upc
export LD_LIBRARY_PATH=../../star-upc/build/:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=../../star-upc/build/:${DYLD_LIBRARY_PATH}
