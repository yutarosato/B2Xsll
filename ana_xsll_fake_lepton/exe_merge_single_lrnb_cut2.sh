#!/bin/bash

./merge_single_lrnb_cut2 $1 $2 $3 $4 $5 << EOF 1>>  log/log_merge_single_lrnb_cut2_${5}_s0${4}.log 2>> error.log
EOF
