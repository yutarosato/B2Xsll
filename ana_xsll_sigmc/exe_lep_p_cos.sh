#!/bin/bash

./draws_2d_lep_p_cos $1 $2 0 << EOF 1>  log/log_2d_lep_p_cos_lep$1_set$2.log 2>> error.log
EOF
