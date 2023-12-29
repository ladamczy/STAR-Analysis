#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -r build/Plots
cd build

g++ -g $DIR/Plots.cxx -o Plots -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I $DIR/../star-upc/include/ -I $DIR/include -L $DIR/../star-upc/build/ -lstar-upc
