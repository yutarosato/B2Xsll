#!/bin/bash

./draws_1d_pid $1 $2 0 << EOF 1>  log/log_pid_pid$1_par$2.log 2>> error.log
EOF
