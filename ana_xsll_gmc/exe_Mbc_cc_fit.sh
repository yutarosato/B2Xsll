#!/bin/bash

./draws_1d_Mbc_cc_fit $1 0 << EOF 1>  log/log_Mbc_cc_fit_s0$1.log 2>> error.log
EOF
