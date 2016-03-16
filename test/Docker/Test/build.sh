#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

docker build -t thewtex/itkultrasound-test $script_dir
