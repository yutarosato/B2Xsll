#!/bin/bash

./draws_1d_Mbc_bb_rest_peak2 $1 $2 0 << EOF 1>  log/log_Mbc_bb_rest_peak2_s0$1.log 2>> error.log
EOF
