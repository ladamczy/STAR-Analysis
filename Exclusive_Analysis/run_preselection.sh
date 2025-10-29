#!/bin/bash
echo "Input file: $1"
echo "Working directory: $(pwd)"
mkdir -p DataAfterPreselection
./build/Preselection "$1"


