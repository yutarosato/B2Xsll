#!/bin/bash

./draws_1d_diagonal_bvtxcl_true_bcs $1 0 << EOF 1> log/log_pdf_bvtxcl_true_bcs_set$1.log 2>> error.log
EOF
