#!/bin/bash

./draws_2d_q2_theta_eff_fraction_kll_p $1 $2 0 << EOF 1>  log/log_2d_q2_theta_eff_fraction_kll_p_lep$1_set$2.log 2>> error.log
EOF
