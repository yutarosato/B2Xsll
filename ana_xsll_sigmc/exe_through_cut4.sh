#!/bin/bash

./through_cut4 $1 $2 << EOF 1>>  log/log_through_cut4.log 2>> error.log
EOF
