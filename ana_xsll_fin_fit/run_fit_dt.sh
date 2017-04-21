#!/bin/csh -f

bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 0;
bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 1;
bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 2;
bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_sub.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_sub.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_sub.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_sub.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_1_6q2.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_1_6q2_sub.sh 15 1;



