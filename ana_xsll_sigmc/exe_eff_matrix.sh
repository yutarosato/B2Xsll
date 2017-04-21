#!/bin/bash

./eff_matrix $1 $2 $3 $4 $5 0 << EOF 1>  log/log_eff_matrix_lep$1_xsid$2_low$3_set$4-$5.log 2>> error.log
EOF
