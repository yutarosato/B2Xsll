#!/bin/bash

./draws_1d_diagonal_M_xs $1 $2 $3 0 << EOF 1> log/log_diag_M_xs_$1_lep$2_set$3.log 2>> error.log
EOF
