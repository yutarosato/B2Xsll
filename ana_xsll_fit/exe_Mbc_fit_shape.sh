#!/bin/bash

./draws_1d_Mbc_fit_shape $1 $2 0 << EOF 1>  log/log_Mbc_fit_shape_$1bin_$2.log 2>> error.log
EOF
