#!/bin/bash

#processing root file used as an argument
root $(dirname $(readlink $0))/'ROOTtoPDF.cxx("'$1'", "'$2'")' <<EOF
.q
EOF