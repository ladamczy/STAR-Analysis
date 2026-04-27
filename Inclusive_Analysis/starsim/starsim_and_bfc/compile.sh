#!/bin/bash

singularity shell --shell /usr/bin/csh -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif <<EOF

starver pro
cons
rm temp_gccflags.c
EOF
