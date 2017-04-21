#!/bin/bash

./through_cut $1 $2 << EOF 1>> log/log_through_cut.log 2>> error.log
EOF
