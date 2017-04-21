#!/bin/csh -f

( ./mklikelihood $1 $2 $3 >> log/log_mklikelihood_$3.log ) >>& error.log
