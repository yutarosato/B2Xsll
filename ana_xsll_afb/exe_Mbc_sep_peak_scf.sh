#!/bin/bash

./draws_1d_Mbc_sep_peak_scf $1 $2 $3 $4 0 << EOF 1>  log/log_1d_Mbc_sep_peak_scf_s0$1_set$2.log 2>> error.log
EOF
