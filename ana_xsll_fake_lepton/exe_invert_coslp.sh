#!/bin/bash

./invert_coslp $1 $2 << EOF 1>>  log/log_invert_coslp.log 2>> error.log
EOF
