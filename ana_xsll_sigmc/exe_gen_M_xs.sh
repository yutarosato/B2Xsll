#!/bin/bash

./draws_1d_gen_M_xs $1 $2 $3 $4 0 << EOF 1>  log/log_gen_M_xs_lep$1_klong$2_set$3.log 2>> error.log
EOF
