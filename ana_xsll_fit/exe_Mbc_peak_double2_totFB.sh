#!/bin/bash

./draws_1d_Mbc_peak_double2_totFB $1 $2 0 << EOF 1>  log/log_Mbc_peak_double2_totFB_s0$1.log 2>> error.log
EOF
