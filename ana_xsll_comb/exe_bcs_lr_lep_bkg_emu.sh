#!/bin/csh -f

mkdir -p $4

( ./bcs_lr_lep_bkg_emu $1 $2 uds     $3 $4 $5 >> log/log_bcs_lr_lep_emu_bkg_uds_s0${2}_${5}.log     ) >>& error.log
( ./bcs_lr_lep_bkg_emu $1 $2 charm   $3 $4 $5 >> log/log_bcs_lr_lep_emu_bkg_charm_s0${2}_${5}.log   ) >>& error.log
( ./bcs_lr_lep_bkg_emu $1 $2 mixed   $3 $4 $5 >> log/log_bcs_lr_lep_emu_bkg_mixed_s0${2}_${5}.log   ) >>& error.log
( ./bcs_lr_lep_bkg_emu $1 $2 charged $3 $4 $5 >> log/log_bcs_lr_lep_emu_bkg_charged_s0${2}_${5}.log ) >>& error.log
