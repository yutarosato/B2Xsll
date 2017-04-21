#!/bin/bash

./draws_1d_Mbc_fit_bin_dt_mod_PRL_plus_single_nu $1 $2 0 << EOF 1>  log/log_Mbc_bin_dt_mod_PRL_plus_single_nu_func$1_$2q2.log 2>> error.log
EOF
