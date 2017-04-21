#!/bin/bash

./draws_2d_q2_theta_eff_cut4 $1 $2 $3 0 << EOF 1>  log/log_2d_q2_theta_eff_cut4_lep$1_xsid$2_set$3.log 2>> error.log
EOF
