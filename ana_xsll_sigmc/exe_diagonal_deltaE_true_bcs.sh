#!/bin/bash

./draws_1d_diagonal_deltaE_true_bcs $1 $2 0 << EOF 1> log/log_pdf_deltaE_true_bcs_lep$1_set$2.log 2>> error.log
EOF
