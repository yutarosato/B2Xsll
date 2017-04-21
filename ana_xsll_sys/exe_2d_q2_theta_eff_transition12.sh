#!/bin/bash

./draws_2d_q2_theta_eff_transition_12 $1 $2 0 << EOF 1>  log/log_2d_q2_theta_eff_transition12_lep$1_set$2.log 2>> error.log
EOF
