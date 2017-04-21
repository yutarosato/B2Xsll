#!/bin/bash

./calib $1 $2 << EOF 1>>  log/log_calib.log 2>> error.log
EOF
