#!/bin/bash

./draws_1d_var $1 0 << EOF 1>  log/log_var_axis$1.log 2>> error.log
EOF
