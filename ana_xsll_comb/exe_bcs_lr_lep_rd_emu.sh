#!/bin/csh -f

mkdir -p $3

( ./bcs_lr_lep_bkg_emu $1 999 rd $2 $3 $4 >> log/log_bcs_lr_lep_emu_rd_${4}.log ) >>& error.log
