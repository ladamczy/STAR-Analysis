#!/bin/bash


#checking if there is an unusual config
while getopts "c:" flag
do
    case $flag in
        c) configur=${OPTARG};;
    esac
    shift
    shift
done

if [[ -n $configur ]]; then
    echo Loaded config no $configur
fi

#processing root file used as an argument
root $(dirname $(readlink $0))/'ROOTtoPDF.cxx("'$1'", "'$2'", "'$configur'")' <<EOF
.q
EOF