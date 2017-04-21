#!/bin/bash

./draws_1d_Mbc $1 $2 $3 $4 0 << EOF 1>  log/log_1d_Mbc_s0$1_set$2.log 2>> error.log
EOF
