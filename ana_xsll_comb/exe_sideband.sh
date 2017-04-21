#!/bin/bash

./draws_1d_sideband $1 $2 $3 $4 0 << EOF 1>  log/log_bgsup_sideband_lep$4_axis$3_s0$1_set$2.log 2>> error.log
EOF
