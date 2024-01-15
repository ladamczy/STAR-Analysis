#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -r build/ThetaStar_K0s
cd build

g++ -g $DIR/ThetaStar_K0.cxx -o ThetaStar_K0s -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I $DIR/../star-upc/include/ -I $DIR/include -L $DIR/../star-upc/build/ -lstar-upc
