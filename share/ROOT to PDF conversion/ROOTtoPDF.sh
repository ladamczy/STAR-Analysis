#!/bin/bash

#processing root file used as an argument
root 'ROOTtoPDF.cxx("'$1'", "'$2'")' <<EOF
.q
EOF