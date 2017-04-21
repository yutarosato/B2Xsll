#!/bin/bash

./draws_1d_Mbc_fit_bin_dt_mod_PRL_sub $1 $2 0 << EOF 1>  log/log_Mbc_bin_dt_mod_PRL_func$1_$2q2_sub.log 2>> error.log
EOF
