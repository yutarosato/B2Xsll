#!/bin/bash

./draws_1d_Mbc_fit_dt $1 0 << EOF 1>  log/log_Mbc_fit_dt_func$1.log 2>> error.log
EOF
