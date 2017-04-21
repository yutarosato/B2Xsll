#!/bin/bash

./draws_1d_entry_gen $1 $2 $3 $4 0 << EOF 1>  log/log_entry_gen_mc$1_exp$2.log 2>> error.log
EOF
