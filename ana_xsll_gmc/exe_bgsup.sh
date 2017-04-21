#!/bin/bash

./draws_1d_bgsup $1 $2 0 << EOF 1>  log/log_bgsup_a$2_s0$1.log 2>> error.log
EOF
