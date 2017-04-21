#!/bin/bash

./draws_1d_bvtxchi2 $1 $2 0 << EOF 1>  log/log_1d_bvtxchi2_s0$1.log 2>> error.log
EOF
