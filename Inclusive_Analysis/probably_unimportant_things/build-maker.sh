#!/bin/bash

SRC="src"
INC="include"
DEST="StRoot/StUPCFilterMaker"

srclist=( StUPCBemcCluster StUPCEvent StUPCFilterMaker StUPCFilterBemcUtil StUPCFilterTrgUtil
StUPCTrack StUPCVertex StUPCTofHit
StRPEvent StUPCRpsTrack StUPCRpsTrackPoint StUPCRpsCluster StUPCFilterRPUtil )

mkdir -p $DEST

i=0
for file in ${srclist[@]}
do
  cp -f $SRC"/"${srclist[$i]}".cxx" $DEST/
  cp -f $INC"/"${srclist[$i]}".h" $DEST/
  ((i++))
done

#TOF calib maker with start time override
cvs co StRoot/StBTofCalibMaker

cons