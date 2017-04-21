#!/bin/bash

./draws_1d_Mbc_emu $1 $2 $3 $4 0 << EOF 1>  log/log_1d_Mbc_emu_lep$3_func$4_s0$1.log 2>> error.log
EOF
