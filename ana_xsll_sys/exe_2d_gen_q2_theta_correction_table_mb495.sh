#!/bin/bash

./draws_2d_gen_q2_theta_correction_table_mb495 $1 $2 $3 $4 0 << EOF 1>  log/log_2d_gen_q2_theta_correction_table_mb495_lep$1_$2_$3_$4.log 2>> error.log
EOF
