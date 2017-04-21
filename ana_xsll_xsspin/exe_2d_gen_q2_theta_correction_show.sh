#!/bin/bash

./draws_2d_gen_q2_theta_correction_show $1 $2 $3 0 << EOF 1>  log/log_2d_gen_q2_theta_correction_show_lep$1_xs$2_set$3.log 2>> error.log
EOF
