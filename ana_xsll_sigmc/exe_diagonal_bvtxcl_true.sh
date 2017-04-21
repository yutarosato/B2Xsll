#!/bin/bash

./draws_1d_diagonal_bvtxcl_true $1 $2 $3 $4 $5 $6 0 << EOF 1> log/log_diag_bvtxcl_true_lep$1_m$2m_set$6.log 2>> error.log
EOF
