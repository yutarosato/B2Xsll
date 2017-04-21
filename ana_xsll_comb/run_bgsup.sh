#!/bin/csh -f

#make;

#set stream  = '0-5'
#set setname = 'A-U'

#set stream  = '3-5'
#set setname = 'A-J'

#set stream  = '0-2'
#set setname = 'K-U'

set stream  = '3-5'
set setname = 'K-U'

#set list_axis = "6 8 9" # 0(evsi) 1(mmiss) 2(deroe) 3(kfbchi) 4(dzll) 5(cos-theta_B) 6(Fmiss) 7(ksfw) 8(Fmiss_qq) 9(Fmiss_bb)
#set list_axis = "0 1 2"
set list_axis = "4 5"

if( 1 == 1 ) then
  foreach axis ($list_axis)
    bsub -q s ./exe_bgsup.sh  $stream $setname $axis  #  [stream][setname][fl_axis]
  end
endif

#bsub -q s ./exe_diagonal_deltaE_true_bcs.sh  $stream $setname 0 # [stream][setname][fl_mode_ll]
#bsub -q s ./exe_diagonal_deltaE_true_bcs.sh  $stream $setname 1
