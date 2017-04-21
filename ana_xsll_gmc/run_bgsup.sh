#!/bin/csh -f

#make;
set list_axis = "3" # 0(evis) 1(mmiss) 2(deroe) 3(kfbchi) 4(dzll) 5(cos-theta_B) 6(Fmiss)
set stream    = '3-5'


if( 1 == 1 ) then
  foreach axis ($list_axis)
    bsub -q s ./exe_bgsup.sh  $stream $axis  # [stream][fl_axis]
  end
endif

#bsub -q s ./exe_diagonal_deltaE_true_bcs.sh $stream 0 # [stream][fl_mode_ll]
#bsub -q s ./exe_diagonal_deltaE_true_bcs.sh $stream 1

