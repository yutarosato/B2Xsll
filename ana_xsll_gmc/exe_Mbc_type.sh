#!/bin/bash

./draws_1d_Mbc_type $1 $2 0 << EOF 1>  log/log_Mbc_type_s0$1.log 2>> error.log
EOF
