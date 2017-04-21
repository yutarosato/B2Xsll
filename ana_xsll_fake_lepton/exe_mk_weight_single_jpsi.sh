#!/bin/bash

./mk_weight_single_jpsi $1 $2 $3 $4 $5 << EOF 1>>  log/log_mk_weight_single_jpsi_${5}_s0${4}.log 2>> error.log
EOF
