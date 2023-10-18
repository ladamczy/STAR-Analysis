#!/bin/bash

#help
Help()
{
   # Display Help
   echo "Script helping with mass drawing histograms"
   echo "Syntax is described in README.md in share/ROOT_to_PDF_conversion folder"
   echo
}

#initialization
mainconfig=0
oneconfig=0
twoconfig=0
#checking if there is an unusual config
while getopts "m:o:t:h" flag
do
    case $flag in
        m) mainconfig=${OPTARG}
        shift
        shift;;
        o) oneconfig=${OPTARG}
        shift
        shift;;
        t) twoconfig=${OPTARG}
        shift
        shift;;
        h) Help
        exit;;
    esac
done

if [[ -n $mainconfig ]]; then
    echo Loaded main config no $mainconfig
fi
if [[ -n $oneconfig ]]; then
    echo Loaded one-dim histogram config no $oneconfig
fi
if [[ -n $twoconfig ]]; then
    echo Loaded two-dim histogram config no $twoconfig
fi

#processing root file used as an argument
root $(dirname $(readlink $0))/'ROOTtoPDF.cxx("'$1'", "'$2'", '$mainconfig', '$oneconfig', '$twoconfig')' <<EOF
.q
EOF