#!/bin/bash

# setting default value so that only one batch is made in default case
length=$1
length="${length:-0}"

for i in $( eval echo {0..$length} );
do
    echo $i
done