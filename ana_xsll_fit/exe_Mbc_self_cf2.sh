#!/bin/bash

./draws_1d_Mbc_self_cf2 $1 0 << EOF 1>  log/log_Mbc_self_cf2_set$1.log 2>> error.log
EOF
