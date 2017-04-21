#!/bin/bash

./draws_1d_M_ll $1 $2 0 << EOF 1>  log/log_M_ll_lep$1_cut$2.log 2>> error.log
EOF
