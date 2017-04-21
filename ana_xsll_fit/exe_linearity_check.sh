#!/bin/bash

./linearity_check $1 0 << EOF 1>  log/log_linearity_check_$1.log 2>> error.log
EOF
