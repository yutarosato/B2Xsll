#!/bin/bash

./draws_1d_Mbc_mode $1 $2 $3 $4 $5 $6 $7 0 << EOF 1>  log/log_Mbc_m$7m_lep$5_q2_$6_s0$1_set$2.log 2>> error.log
EOF
