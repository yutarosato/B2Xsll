#!/bin/bash

./mk_weight_single_jpsi_data $1 $2 $3 << EOF 1>>  log/log_mk_weight_single_jpsi_data.log 2>> error.log
EOF
