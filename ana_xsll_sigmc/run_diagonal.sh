#!/bin/bash

#make;

setname="K-U" # A, AB, A-C


###############################################################
# plot distributions
###############################################################
if [ 1 == 0 ]; then
  for xs in 1 10  101 1001 110 1010  201 1101 210 1110  301 1201 310 1210  401 1301 410 1310
  do
    #[mu]
      bsub -q s ./exe_diagonal_Mbc.sh       $xs 0 $setname
      bsub -q s ./exe_diagonal_deltaE.sh    $xs 0 $setname
      bsub -q s ./exe_diagonal_M_xs.sh      $xs 0 $setname
      bsub -q s ./exe_diagonal_M_ll.sh      $xs 0 $setname
    #[e]
      bsub -q s ./exe_diagonal_Mbc.sh       $xs 1 $setname
      bsub -q s ./exe_diagonal_deltaE.sh    $xs 1 $setname
      bsub -q s ./exe_diagonal_M_xs.sh      $xs 1 $setname
      bsub -q s ./exe_diagonal_M_ll.sh      $xs 1 $setname
  done
fi

      
###############################################################
# check mode dependence of distribution
###############################################################
#bsub -q s ./exe_diagonal_deltaE_true.sh  1   1    10  101  110 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  1 201   210  101  110 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  1 301   310  101  110 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  1 401   410  101  110 $setname

#bsub -q s ./exe_diagonal_deltaE_true.sh  1 1101 1110 1001 1010 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  1 1201 1210 1001 1010 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  1 1301 1310 1001 1010 $setname

#bsub -q s ./exe_diagonal_deltaE_true.sh  0    1   10  101  110 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  0  201  210  101  110 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  0  301  310  101  110 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  0  401  410  101  110 $setname

#bsub -q s ./exe_diagonal_deltaE_true.sh  0 1101 1110 1001 1010 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  0 1201 1210 1001 1010 $setname
#bsub -q s ./exe_diagonal_deltaE_true.sh  0 1301 1310 1001 1010 $setname

#bsub -q s ./exe_diagonal_bvtxcl_true.sh  0 1301 1310 1001 1010 $setname
      
###############################################################
# PDF determination
###############################################################
#bsub -q s ./exe_diagonal_deltaE_true_bcs.sh 0 $setname
#bsub -q s ./exe_diagonal_deltaE_true_bcs.sh 1 $setname
bsub -q s ./exe_diagonal_bvtxcl_true_bcs.sh   $setname
