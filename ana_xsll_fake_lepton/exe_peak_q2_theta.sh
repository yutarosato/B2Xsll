#!/bin/bash

./draws_2d_peak_q2_theta $1 0 << EOF 1>  log/log_peak$1_q2_theta.log 2>> error.log
EOF
