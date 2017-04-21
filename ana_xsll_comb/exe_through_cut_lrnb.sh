#!/bin/bash

./through_cut_lrnb $1 $2 << EOF 1>>  log/log_through_cut_lrnb.log 2>> error.log
EOF
