#!/bin/bash

./draws_1d_all $1 $2 $3 $4 0 << EOF 1>  log/log_all_$1_lep$2_gmxs$3_set$4.log 2>> error.log
EOF
