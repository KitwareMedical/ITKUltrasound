#!/bin/bash

# This is a script to build the modules and run the test suite in the base
# Docker container.

die() {
  echo "Error: $@" 1>&2
  exit 1;
}

cd /usr/src/ITKUltrasound-build || die "Could not cd into the build directory"

cmake \
  -G Ninja \
  -DITK_DIR:PATH=/usr/src/ITK-build \
  -DITKUltrasound_USE_VTK:BOOL=ON \
  -DPYTHON_EXECUTABLE:FILEPATH=/usr/bin/python3 \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILDNAME:STRING=External-Ultrasound \
    /usr/src/ITKUltrasound || die "CMake configuration failed"
ctest -VV -D Experimental || die "ctest failed"
