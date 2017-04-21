#!/bin/csh -f

(./bcs_merge $1 $2 $3 $4 >> log/log_bcs_merge_set$2.log ) >>& error.log
