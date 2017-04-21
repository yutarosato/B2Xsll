#!/bin/bash

./draws_1d_ksfw_var_mode $1 $2 $3 $4 0 << EOF 1>  log/log_ksfw_var_lep$3_ksfw_$4_s0$1_set$2.log 2>> error.log
EOF
