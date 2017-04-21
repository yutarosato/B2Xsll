#!/bin/bash

./draws_1d_diagonal_deltaE_true_bcs $1 $2 $3 0 << EOF 1>  log/log_diagonal_deltaE_true_bcs_s0$1_set$2_lep$3.log 2>> error.log
EOF
