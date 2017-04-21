#!/bin/bash

./draws_2d_q2_theta_gmc1 $1 $2 $3 0 << EOF 1>  log/log_2d_q2_theta_gmc1_lep$3_s0$1.log 2>> error.log
EOF
