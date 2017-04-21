#!/bin/bash

./cal_delh $1 $2 << EOF 1>>  log/log_cal_delh.log 2>> error.log
EOF
