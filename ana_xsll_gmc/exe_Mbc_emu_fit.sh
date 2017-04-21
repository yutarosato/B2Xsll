#!/bin/bash

./draws_1d_Mbc_emu_fit $1 $2 0 << EOF 1>  log/log_Mbc_emu_fit_func$2_s0$1.log 2>> error.log
EOF
