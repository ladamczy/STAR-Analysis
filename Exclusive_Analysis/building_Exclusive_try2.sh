#!/bin/bash

# Check ROOT version
REQUIRED_ROOT_VERSION="6.34"
CURRENT_ROOT_VERSION=$(root-config --version | cut -d. -f1-2)

echo $ROOTSYS
root-config --prefix

if [[ "$CURRENT_ROOT_VERSION" != "$REQUIRED_ROOT_VERSION" ]]; then
    echo "Error: Required ROOT version is $REQUIRED_ROOT_VERSION, but found $CURRENT_ROOT_VERSION"
    exit 1
fi

# Check if build directory exists
if [ -d "build" ]; then
    echo "Cleaning build directory..."
    rm -rf build/*
else
    echo "Creating build directory..."
    mkdir -p build
fi

# Get ROOT flags and libraries
ROOT_FLAGS=$(root-config --cflags)
ROOT_LIBS=$(root-config --libs)
ROOT_EG_LIB="-lEG"  # Add EG library explicitly

# Set the library path
export LD_LIBRARY_PATH=/home/sbhosale/Work/STAR-Analysis/star-upc-new/build:$ROOTSYS/lib

# Compile Exclusive_try2
echo "Compiling Exclusive_try2..."
    g++ Exclusive_try2.cxx \
        ExclusiveCode.cxx \
        -I/usr/include/root \
        -I./include \
        -I/home/sbhosale/Work/STAR-Analysis/star-upc-new/include \
        -pthread -std=c++17 -m64 \
        $(root-config --libs) -lEG \
        -L/home/sbhosale/Work/STAR-Analysis/star-upc-new/build -lstar-upc


if [ $? -eq 0 ]; then
    echo "Build complete. Executable is in build/Preselection"
else
    echo "Build failed!"
    exit 1
fi
