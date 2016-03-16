#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

docker build -t thewtex/itkultrasound-base $script_dir
