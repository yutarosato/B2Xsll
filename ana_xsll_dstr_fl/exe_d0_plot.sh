#!/bin/bash

./draws_1d_d0_plot $1 $2 $3 $4 $5 0 << EOF 1>  log/log_d0_plot_pid$1_par$2_mom$3_cos$4_adef$5.log 2>> error.log
EOF
