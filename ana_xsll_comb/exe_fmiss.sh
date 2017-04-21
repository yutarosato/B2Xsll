#!/bin/bash

./draws_1d_fmiss $1 $2 $3 0 << EOF 1>  log/log_fmiss_axis$3_s0$1_set$2.log 2>> error.log
EOF
