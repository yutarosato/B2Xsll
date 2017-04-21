#!/bin/bash

./linearity_check_bin $1 0 << EOF 1>  log/log_linearity_check_bin_$1.log 2>> error.log
EOF
