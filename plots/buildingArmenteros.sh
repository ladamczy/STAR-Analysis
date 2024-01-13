#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -r build/Podolanski-Armenteros
cd build

g++ -g $DIR/Podolanski-Armenteros.cxx -o Podolanski-Armenteros -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I $DIR/../star-upc/include/ -I $DIR/include -L $DIR/../star-upc/build/ -lstar-upc
