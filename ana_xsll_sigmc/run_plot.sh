#!/bin/csh -f

#make;

#set setname = 'A'   # small dataset for debug
set setname = 'A-J' # official dataset for significance

if( 0 == 1 ) then
  foreach axis (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)
    bsub -q e ./exe_theta_ksfw.sh  $setname $axis  #  [setname][fl_axis]
  end
endif

if( 1 == 1 ) then
  foreach xs (1 10    101 1001 110 1010    201 1101 210 1110    301 1201 310 1210    401 1301 410 1310)
    bsub -q e ./exe_ksfw.sh  $xs 1 $setname  #  [fl_mode_xs][fl_mode_ll][setname]
    bsub -q e ./exe_ksfw.sh  $xs 0 $setname  #  [fl_mode_xs][fl_mode_ll][setname]
  end
endif
