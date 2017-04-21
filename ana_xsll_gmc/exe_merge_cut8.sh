#!/bin/bash

mkdir -p $2

./merge_cut8 $1 $2 $3 $4 $5 << EOF 1>>  log/log_merge_cut8_${5}_s0${4}.log 2>> error.log
EOF
