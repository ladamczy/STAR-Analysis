#!/usr/bin/bash

# cd ~/star-upc
cd $(dirname $0)
# cp ./Different_CMakes/CMakeListsNormal.txt ./CMakeLists.txt
# rm -r build
# mkdir build
cd build
cmake ..
make