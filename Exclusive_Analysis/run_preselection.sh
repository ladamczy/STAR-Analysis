#!/bin/bash

echo "Input file: $1"
echo "Working directory: $(pwd)"
echo "Hostname: $(hostname)"

# Create output directory with absolute path
OUTPUT_DIR="/home/sbhosale/Work/STAR-Analysis/Exclusive_Analysis/DataAfterPreselection"
mkdir -p "$OUTPUT_DIR"

# Run the program
./build/Preselection "$1"

# List what was created
echo "Output files created:"
ls -lh "$OUTPUT_DIR"
ls -l 
