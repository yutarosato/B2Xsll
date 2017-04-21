#!/bin/bash

./through_cut3 $1 $2 << EOF 1>>  log/log_through_cut3.log 2>> error.log
EOF
