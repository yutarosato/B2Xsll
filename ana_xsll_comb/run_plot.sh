#!/bin/csh -f

#make;

#set stream  = '0'   # small dataset for debug
#set nstream = '1'   # small dataset for debug
#set setname = 'A'   # small dataset for debug
#set nset    = '1'   # small dataset for debug
set stream  = '0-2' # official dataset for significance
set nstream = '3'   # official dataset for significance
set setname = 'A-J' # official dataset for significance
set nset    = '10'  # official dataset for significance
     
###########################################################################################3
if( 1 == 0 ) then
  foreach axis (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)
    bsub -q s ./exe_ksfw.sh  $stream $setname $axis  #  [stream][setname][fl_axis]
  end
endif

  
###########################################################################################3
if( 1 == 0 ) then
  foreach axis (0 1 2 3 4 5)
    bsub -q s ./exe_fmiss.sh  $stream $setname $axis  #  [stream][setname][fl_axis]
  end
endif

  
###########################################################################################3
#bsub -q s ./exe_ksfw_lr.sh  $stream $setname 1  0 # [stream][setname][fl_mode_ll][fl_xs_region] 
#bsub -q s ./exe_ksfw_lr.sh  $stream $setname 1 -1 
#bsub -q s ./exe_ksfw_lr.sh  $stream $setname 1  1



if( 0 == 1 ) then
  foreach ksfw (0 1 2 3 4 5 6)
     bsub -q s ./exe_ksfw_var.sh  $stream $setname 1   0 $ksfw  # [stream][setname][fl_mode_ll][fl_xs_region][fl_ksfw_region]
     bsub -q s ./exe_ksfw_var.sh  $stream $setname 1  -1 $ksfw
     bsub -q s ./exe_ksfw_var.sh  $stream $setname 1   1 $ksfw
     bsub -q s ./exe_ksfw_var.sh  $stream $setname 0   0 $ksfw
     bsub -q s ./exe_ksfw_var.sh  $stream $setname 0  -1 $ksfw
     bsub -q s ./exe_ksfw_var.sh  $stream $setname 0   1 $ksfw
  end
endif

if( 1 == 1 ) then
  foreach ksfw (0 1 2 3 4 5 6)
     bsub -q s ./exe_ksfw_var_mode.sh  $stream $setname 1 $ksfw  # [stream][setname][fl_mode_ll][fl_xs_region][fl_ksfw_region]
     bsub -q s ./exe_ksfw_var_mode.sh  $stream $setname 0 $ksfw
  end
endif
  
###########################################################################################3

if( 0 == 1 ) then
  foreach axis (0 1 2 3 4 5 6 7) # 0(evis1), 1(abs(mmiss)), 2(deroe), 3(kfbchi/kdbdgf), 4(dzll3d), 5(bccm), 6(fmiss) 7(ksfw)
     bsub -q s ./exe_sideband.sh  $stream $setname $axis 0  # [stream][setname][fl_axis][fl_mode_ll]
     bsub -q s ./exe_sideband.sh  $stream $setname $axis 1
  end
endif
  

