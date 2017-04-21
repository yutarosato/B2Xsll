#!/bin/bash

./draws_1d_Mbc_mode $1 $2 $3 0 << EOF 1>  log/log_Mbc_mode_m$1m_mc$2_set$3.log 2>> error.log
EOF
