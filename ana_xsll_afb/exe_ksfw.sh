#!/bin/bash

./draws_1d_ksfw $1 $2 $3 0 << EOF 1>  log/log_1d_ksfw_axis$1_s0$2.log 2>> error.log
EOF
