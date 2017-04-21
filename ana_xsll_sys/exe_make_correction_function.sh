#!/bin/bash

./make_correction_function $1 0 << EOF 1>  log/log_make_correction_function_axis$1.log 2>> error.log
EOF
