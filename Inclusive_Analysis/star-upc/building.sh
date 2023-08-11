#!/usr/bin/bash

# cd ~/star-upc
cd $(dirname $0)
rm -r build
mkdir build
cd build
cmake ..
make