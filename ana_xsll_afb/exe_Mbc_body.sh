#!/bin/bash

./draws_1d_Mbc_body $1 $2 $3 $4 $5 0 << EOF 1>  log/log_1d_Mbc_$5body_s0$1_set$2.log 2>> error.log
EOF
