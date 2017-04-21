#!/bin/csh -f

#make;

#set setname = 'A'   # small dataset for debug
set setname = 'A-J' # official dataset for significance

bsub -q e ./exe_lep_p_cos.sh  0 $setname  # [fl_mode_ll] [setname]
bsub -q e ./exe_lep_p_cos.sh  1 $setname  #
