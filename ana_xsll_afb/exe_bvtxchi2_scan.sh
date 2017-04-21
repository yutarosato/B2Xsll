#!/bin/bash

./draws_1d_bvtxchi2_scan $1 $2 $3 0 << EOF 1>  log/log_1d_bvtxchi2_scan_factor$1_s0$2.log 2>> error.log
EOF
