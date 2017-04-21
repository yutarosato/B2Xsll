#!/bin/bash

./draws_1d_Mbc_nopeak_fit $1 $2 $3 0 << EOF 1>  log/log_Mbc_nopeak_fit_func$3_s0$1.log 2>> error.log
EOF
