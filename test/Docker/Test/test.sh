#!/bin/bash

# This is a script to build the modules and run the test suite in the base
# Docker container.

set -x
set -o

cd /usr/src/ITKUltrasound
branch=$(git rev-parse --abbrev-ref HEAD)
date=$(date +%F_%H_%M_%S)
sha=$(git rev-parse --short HEAD)

cd /usr/src/ITKUltrasound-build

cmake \
  -G Ninja \
  -DITK_DIR:PATH=/usr/src/ITK-build \
  -DITKUltrasound_USE_VTK:BOOL=ON \
  -DPYTHON_EXECUTABLE:FILEPATH=/opt/conda/bin/python \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILDNAME:STRING=External-Ultrasound-${branch}-${date}-${sha} \
    /usr/src/ITKUltrasound
ctest -VV -D Experimental
