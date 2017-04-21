#!/bin/bash

./draws_1d_entry $1 $2 $3 0 << EOF 1> log/log_entry_lep$1_gmxs$2_set$3.log 2>> error.log
EOF
