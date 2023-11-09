#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -r build/MonteCarloEfficiency
cd build

#export LD_LIBRARY_PATH=../../star-upc/build/:$ROOTSYS/lib
g++  ../MonteCarloEfficiency.cxx  -o MonteCarloEfficiency -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I ../../star-upc/include/ -I ../include -L ../../star-upc/build/ -lstar-upc
