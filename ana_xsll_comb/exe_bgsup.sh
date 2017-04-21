#!/bin/bash

./draws_1d_bgsup $1 $2 $3 0 << EOF 1>  log/log_bgsup_a$3_s0$1_set$2.log 2>> error.log
EOF
