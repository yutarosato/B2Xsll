#!/bin/bash

./through_single_cut2 $1 $2 << EOF 1>>  log/log_through_single_cut2.log 2>> error.log
EOF
