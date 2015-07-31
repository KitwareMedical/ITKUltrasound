#!/bin/sh

# The default password is "wavy".  To use a different password, pass it in
# as an argument to this script.
password_arg=""
if [ $# -gt 1 ]; then
  password_arg="-e PASSWORD=$1"
fi

script_dir="`cd $(dirname $0); pwd`"

docker run \
  --rm \
  -p 443:8888 \
  $password_arg \
  -v $script_dir/../../examples/Notebooks/:/notebooks/:rw \
  itkultrasound-base
