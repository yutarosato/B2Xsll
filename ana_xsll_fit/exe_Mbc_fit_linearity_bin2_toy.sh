#!/bin/bash

./draws_1d_Mbc_fit_linearity_bin2_toy $1 $2 $3 $4 $5 $6 $7 $8 0 << EOF 1>  log/log_Mbc_bin2_toy_func$5_ratio$6_s0$1_set$2_$7q2_$8.log 2>> error.log
EOF
