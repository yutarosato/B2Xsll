#!/bin/bash

./draws_2d_theta_ksfw $1 $2 0 << EOF 1>  log/log_theta_ksfw_axis$2_set$1.log 2>> error.log
EOF
