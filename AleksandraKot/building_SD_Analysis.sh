#!/bin/bash

# Check if build directory exists
if [ -d "build" ]; then
    echo "Cleaning build directory..."
    rm -rf build/SD_Analysis
else
    echo "Creating build directory..."
    mkdir -p build
fi

# Get ROOT flags and libraries
ROOT_FLAGS=$(root-config --cflags)
ROOT_LIBS=$(root-config --libs)
ROOT_EG_LIB="-lEG"  # Add EG library explicitly

# Define the library directory
STAR_UPC_LIB_DIR="../star-upc-new/build"

# Compile Preselection
echo "Compiling SD_Analysis..."

g++ SD_Analysis.cxx \
    -o build/SD_Analysis \
    -I/usr/include/root \
    -I./include \
    -I${STAR_UPC_LIB_DIR%/build}/include \
    -pthread -std=c++17 -m64 \
    ${ROOT_FLAGS} ${ROOT_LIBS} ${ROOT_EG_LIB} \
    -L${STAR_UPC_LIB_DIR} -lstar-upc \
    -Wl,-rpath,${STAR_UPC_LIB_DIR}

if [ $? -eq 0 ]; then
    echo "Build complete. Executable is in build/SD_Analysis"
else
    echo "Build failed!"
    exit 1
fi
