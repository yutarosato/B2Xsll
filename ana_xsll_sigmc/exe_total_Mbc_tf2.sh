#!/bin/bash

./draws_1d_total_Mbc_tf2 $1 $2 0 << EOF 1> log/log_total_Mbc_tf2_lep$1_set$2.log 2>> error.log
EOF
