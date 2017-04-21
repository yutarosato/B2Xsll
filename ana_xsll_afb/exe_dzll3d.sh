#!/bin/bash

./draws_1d_dzll3d $1 $2 0 << EOF 1>  log/log_1d_dzll3d_s0$1.log 2>> error.log
EOF
