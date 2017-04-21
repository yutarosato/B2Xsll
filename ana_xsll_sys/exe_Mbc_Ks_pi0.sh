#!/bin/bash

./draws_1d_Mbc_Ks_pi0 $1 $2 $3 0 << EOF 1>  log/log_Mbc_Ks_pi0_m$1m_mc$2_set$3.log 2>> error.log
EOF
