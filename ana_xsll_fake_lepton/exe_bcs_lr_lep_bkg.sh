#!/bin/csh -f

mkdir -p $2

( ./bcs_lr_lep_bkg $1 $2 $3 >> log/log_bcs_lr_lep_bkg.log ) >>& error.log
