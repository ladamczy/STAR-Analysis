#!/usr/bin/bash

cd ~/star-upc
rm -r build
mkdir build
cd build
cmake ..
make