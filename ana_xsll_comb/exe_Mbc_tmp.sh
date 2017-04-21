#!/bin/bash

./draws_1d_Mbc_tmp $1 $2 $3 $4 $5 0 << EOF 1>  log/log_Mbc_q2_$5_s0$1_set$2.log 2>> error.log
EOF
