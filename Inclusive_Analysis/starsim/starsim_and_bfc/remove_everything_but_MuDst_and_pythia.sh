#!/bin/bash

if [ "$1" -eq "1" ]; then
    echo "Step 1"
    ls -lh

    rename .MuDst.root .save *.*
    rename .root .delete *.*.root

    echo "Step 2"
    ls -lh

    rm *.delete
    rename .save .MuDst.root *.*

    echo "Result"
    ls -lh
fi
