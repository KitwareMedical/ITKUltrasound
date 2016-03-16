#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

docker run \
  -v $script_dir/../../..:/usr/src/ITKUltrasound \
  thewtex/itkultrasound-test \
    /usr/src/ITKUltrasound/test/Docker/Test/test.sh
