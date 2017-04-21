#!/bin/bash

#make;

#[mu]
#bsub -q e ./exe_eff_matrix.sh 0 0 1  0  9; # [fl_mode_ll][fl_xsid][fl_low][setname_i][setname_f]
#bsub -q e ./exe_eff_matrix.sh 0 0 1 10 20;
#bsub -q e ./exe_eff_matrix.sh 0 0 1  0 20;
#bsub -q e ./exe_eff_matrix.sh 0 1 1  0 20;
#bsub -q e ./exe_eff_matrix.sh 0 2 1  0 20;
#bsub -q e ./exe_eff_matrix.sh 0 3 1  0 20;
bsub -q e ./exe_eff_matrix.sh 0 0 1  1  1;

#[e]
#bsub -q e ./exe_eff_matrix.sh 1 0 1  0  9;
#bsub -q e ./exe_eff_matrix.sh 1 0 1 10 20;
#bsub -q e ./exe_eff_matrix.sh 1 0 1  0 20;
#bsub -q e ./exe_eff_matrix.sh 1 1 1  0 20;
#bsub -q e ./exe_eff_matrix.sh 1 2 1  0 20;
#bsub -q e ./exe_eff_matrix.sh 1 3 1  0 20;

bsub -q e ./exe_eff_matrix.sh 1 0 1  1  1;
