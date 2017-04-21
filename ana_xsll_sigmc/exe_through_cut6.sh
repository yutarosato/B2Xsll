#!/bin/bash

./through_cut6 $1 $2 << EOF 1>>  log/log_through_cut6.log 2>> error.log
EOF
