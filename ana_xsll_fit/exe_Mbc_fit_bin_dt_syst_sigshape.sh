#!/bin/bash

./draws_1d_Mbc_fit_bin_dt_syst_sigshape $1 $2 $3 0 << EOF 1>  log/log_Mbc_bin_dt_syst_sigshape_func$1_$2q2_cnt$3.log 2>> error.log
EOF
