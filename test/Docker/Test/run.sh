#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

docker run \
  --rm \
  -v $script_dir/../../..:/usr/src/ITKUltrasound \
  kitwaremedical/itkultrasound-test \
    /usr/src/ITKUltrasound/test/Docker/Test/test.sh
