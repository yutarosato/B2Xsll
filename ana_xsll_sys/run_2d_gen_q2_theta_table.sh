#!/bin/bash

#bsub -q e ./exe_2d_gen_q2_theta_correction_table.sh                    1  norm  100 100
#bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kll_p.sh     1  norm  100 100
#bsub -q e ./exe_2d_gen_q2_theta_correction_table.sh                    0  norm  100 100
#bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kll_p.sh     0  norm  100 100


# [TABLE]
#for A7  in norm # norm flip
#do
#for A9  in `seq -200 20 200`
#do
#for A10 in `seq -200 20 200`
#do
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table.sh                    1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table.sh                    0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fm200.sh              1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fm200.sh              0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fm480.sh              1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fm480.sh              0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_mb465.sh              1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_mb465.sh              0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_mb495.sh              1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_mb495.sh              0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_hadronization.sh      1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_hadronization.sh      0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_transition10.sh       1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_transition10.sh       0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_transition12.sh       1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_transition12.sh       0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kll_p.sh     1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kll_p.sh     0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kll_m.sh     1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kll_m.sh     0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kstrll_p.sh  1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kstrll_p.sh  0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kstrll_m.sh  1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_kstrll_m.sh  0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_xsll_p.sh    1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_xsll_p.sh    0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_xsll_m.sh    1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_fraction_xsll_m.sh    0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_lambdaone429.sh        1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_lambdaone429.sh        0  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_lambdaone362.sh        1  $A7  $A9  $A10
   #bsub -q e ./exe_2d_gen_q2_theta_correction_table_lambdaone362.sh        0  $A7  $A9  $A10
#done
#done
#done

# [MAKE TEXT DATA]
#grep -h HOGE correction_table_6qbin_with_corr/log_*                    > correction_table.dat;                     echo "*" >> correction_table.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fm200/log_*              > correction_table_fm200.dat;               echo "*" >> correction_table_fm200.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fm480/log_*              > correction_table_fm480.dat;               echo "*" >> correction_table_fm480.dat;
#grep -h HOGE correction_table_6qbin_with_corr_mb465/log_*              > correction_table_mb465.dat;               echo "*" >> correction_table_mb465.dat;
#grep -h HOGE correction_table_6qbin_with_corr_mb495/log_*              > correction_table_mb495.dat;               echo "*" >> correction_table_mb495.dat;
#grep -h HOGE correction_table_6qbin_with_corr_lambdaone429/log_*        > correction_table_lambdaone429.dat;        echo "*" >> correction_table_lambdaone429.dat;
#grep -h HOGE correction_table_6qbin_with_corr_lambdaone362/log_*        > correction_table_lambdaone362.dat;        echo "*" >> correction_table_lambdaone362.dat;
#grep -h HOGE correction_table_6qbin_with_corr_transition10/log_*       > correction_table_transition10.dat;        echo "*" >> correction_table_transition10.dat;
#grep -h HOGE correction_table_6qbin_with_corr_transition12/log_*       > correction_table_transition12.dat;        echo "*" >> correction_table_transition12.dat;
#grep -h HOGE correction_table_6qbin_with_corr_hadronization/log_*      > correction_table_hadronization.dat;       echo "*" >> correction_table_hadronization.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fraction_kll_p/log_*     > correction_table_fraction_kll_p.dat;      echo "*" >> correction_table_fraction_kll_p.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fraction_kll_m/log_*     > correction_table_fraction_kll_m.dat;      echo "*" >> correction_table_fraction_kll_m.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fraction_kstrll_p/log_*  > correction_table_fraction_kstrll_p.dat;   echo "*" >> correction_table_fraction_kstrll_p.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fraction_kstrll_m/log_*  > correction_table_fraction_kstrll_m.dat;   echo "*" >> correction_table_fraction_kstrll_m.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fraction_xsll_p/log_*    > correction_table_fraction_xsll_p.dat;     echo "*" >> correction_table_fraction_xsll_p.dat;
#grep -h HOGE correction_table_6qbin_with_corr_fraction_xsll_m/log_*    > correction_table_fraction_xsll_m.dat;     echo "*" >> correction_table_fraction_xsll_m.dat;

#grep -h HOGE correction_table_6qbin_with_corr_test_kll/log_*            > correction_table_test_kll.dat;            echo "*" >> correction_table_test_kll.dat;
#grep -h HOGE correction_table_6qbin_with_corr_test_kstrll/log_*         > correction_table_test_kstrll.dat;         echo "*" >> correction_table_test_kstrll.dat;
#grep -h HOGE correction_table_6qbin_with_corr_test_xsll/log_*           > correction_table_test_xsll.dat;           echo "*" >> correction_table_test_xsll.dat;


#for axis  in `seq 0 1 13`
#for axis  in `seq 14 1 16`
for axis  in 17
do
bsub -q e ./exe_make_correction_function.sh $axis ;
done
