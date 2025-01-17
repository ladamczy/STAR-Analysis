#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd ../star-upc
./buildLib.sh
cd $DIR

rm -r build
mkdir build
cd build
export LD_LIBRARY_PATH=../../star-upc/build/:$ROOTSYS/lib
g++ /home/sbhosale/Work/STAR-Analysis/Exclusive_Analysis/src/MatchFillPosition.cxx /home/sbhosale/Work/STAR-Analysis/Exclusive_Analysis/src/ReadFillPositionFile.cxx ../ExclusiveAnalysisStUPCV0.cxx -o ExclusiveAnalysisStUPCV0 -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I /home/sbhosale/Work/STAR-Analysis/star-upc/include/ -I /home/sbhosale/Work/STAR-Analysis/Exclusive_Analysis/include -L /home/sbhosale/Work/STAR-Analysis/star-upc/build/ -lstar-upc
#g++ ../src/MatchFillPosition.cxx ../include/* ../src/ReadFillPositionFile.cxx ../$1.cxx -o $1 -Wl,--copy-dt-needed-entries `root-config --cflags` `root-config --libs` -I ../../star-upc/include/ -I ../include -L ../../star-upc/build/ -lstar-upc
