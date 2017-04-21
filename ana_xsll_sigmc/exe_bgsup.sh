#!/bin/bash

./draws_1d_bgsup $1 $2 0 << EOF 1>  log/log_pdf_bgsup_$1_set$2.log 2>> error.log
EOF
