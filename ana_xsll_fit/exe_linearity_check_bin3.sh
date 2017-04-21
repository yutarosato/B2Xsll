#!/bin/bash

./linearity_check_bin3 $1 0 << EOF 1>  log/log_linearity_check_bin3_$1.log 2>> error.log
EOF
