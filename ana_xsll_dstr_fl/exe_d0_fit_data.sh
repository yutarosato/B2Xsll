#!/bin/bash

./draws_1d_d0_fit_data $1 $2 $3 $4 0 << EOF 1>  log/log_d0_fit_data_par$1_mom$2_cos$3_adef$4.log 2>> error.log
EOF
