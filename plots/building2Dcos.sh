#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -r build/cos2DK0s
cd build

g++ -g $DIR/cos2DK0s.cxx -o cos2DK0s -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I $DIR/../star-upc/include/ -I $DIR/include -L $DIR/../star-upc/build/ -lstar-upc
