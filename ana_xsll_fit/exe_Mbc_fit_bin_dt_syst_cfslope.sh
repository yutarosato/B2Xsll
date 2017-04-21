#!/bin/bash

./draws_1d_Mbc_fit_bin_dt_syst_cfslope $1 $2 $3 0 << EOF 1>  log/log_Mbc_bin_dt_syst_cfslope_func$1_axis$3_$2q2.log 2>> error.log
EOF
