#!/bin/bash

./through_cut5 $1 $2 << EOF 1>>  log/log_through_cut5.log 2>> error.log
EOF
