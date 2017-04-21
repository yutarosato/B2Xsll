#!/bin/bash

./draws_1d_Mbc_fit_bin_dt_totFB $1 0 << EOF 1>  log/log_Mbc_bin_dt_totFB_func$1.log 2>> error.log
EOF
