#!/bin/bash

./draws_2d_p_cos $1 $2 $3 0 << EOF 1>  log/log_p_cos_pid$1_par$2_adef$3.log 2>> error.log
EOF
