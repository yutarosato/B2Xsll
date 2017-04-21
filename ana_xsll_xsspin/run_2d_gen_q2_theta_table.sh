#!/bin/bash

# [SHOW]
#bsub -q e ./exe_2d_gen_q2_theta_correction_show.sh  0 3 A
#bsub -q e ./exe_2d_gen_q2_theta_correction_show.sh  1 3 A

# [TABLE]
#for A7  in norm flip # norm flip
#do
#for A9  in `seq -200 20 200`
#do
#for A10 in `seq -200 20 200`
#do
#   bsub -q e ./exe_2d_gen_q2_theta_correction_table.sh  1  $A7  $A9  $A10
#   bsub -q e ./exe_2d_gen_q2_theta_correction_table.sh  0  $A7  $A9  $A10
#    echo "./exe_2d_gen_q2_theta_correction_table.sh  1  $A7  $A9  $A10" >> tmp.list
#    echo "./exe_2d_gen_q2_theta_correction_table.sh  0  $A7  $A9  $A10" >> tmp.list
#done
#done
#done

# [MAKE TEXT DATA]
#grep -h HOGE correction_table_6qbin_with_corr_xsspin/log_*         > correction_table_6qbin_with_corr_xsspin.dat;         echo "*" >> correction_table_6qbin_with_corr_xsspin.dat;
#grep -h HOGE correction_table_6qbin_with_corr_xsspin_fr10/log_*     > correction_table_6qbin_with_corr_xsspin_fr10.dat;    echo "*" >> correction_table_6qbin_with_corr_xsspin_fr10.dat;
#grep -h HOGE correction_table_6qbin_with_corr_xsspin_fr50/log_*     > correction_table_6qbin_with_corr_xsspin_fr50.dat;    echo "*" >> correction_table_6qbin_with_corr_xsspin_fr50.dat;
#grep -h HOGE correction_table_6qbin_with_corr_xsspin_fr-10/log_*    > correction_table_6qbin_with_corr_xsspin_fr-10.dat;   echo "*" >> correction_table_6qbin_with_corr_xsspin_fr-10.dat;
#grep -h HOGE correction_table_6qbin_with_corr_xsspin_fr-50/log_*    > correction_table_6qbin_with_corr_xsspin_fr-50.dat;   echo "*" >> correction_table_6qbin_with_corr_xsspin_fr-50.dat;
grep -h HOGE correction_table_6qbin_with_corr_xsspin_fr100/log_*     > correction_table_6qbin_with_corr_xsspin_fr100.dat;    echo "*" >> correction_table_6qbin_with_corr_xsspin_fr100.dat;
grep -h HOGE correction_table_6qbin_with_corr_xsspin_fr-100/log_*    > correction_table_6qbin_with_corr_xsspin_fr-100.dat;   echo "*" >> correction_table_6qbin_with_corr_xsspin_fr-100.dat;
