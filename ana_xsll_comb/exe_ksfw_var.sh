#!/bin/bash

./draws_1d_ksfw_var $1 $2 $3 $4 $5 0 << EOF 1>  log/log_ksfw_var_lep$3_xs_$4_ksfw_$5_s0$1_set$2.log 2>> error.log
EOF
