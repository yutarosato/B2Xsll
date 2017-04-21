#!/bin/bash

./through_cut_lrnb_emu $1 $2 $3 << EOF 1>>  log/log_through_cut_lrnb_emu$3.log 2>> error.log
EOF
