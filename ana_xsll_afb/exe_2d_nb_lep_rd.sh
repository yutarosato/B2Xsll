#!/bin/bash

./draws_2d_nb_lep_rd $1 $2 $3 0 << EOF 1>  log/log_2d_nb_lep$3_s0$1.log 2>> error.log
EOF
