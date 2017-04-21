#!/bin/csh -f

#make;

#set stream  = '0'   # test
#set nstream = '1'   # test
     
#set stream  = '3-5' # pdf
#set nstream = '3'   # pdf
     
#set stream  = '0-2' # evaluation
#set nstream = '3'   # evaluation
     
set stream  = '0-5' # all data
set nstream = '6'   # all data

#*************************************************************************************

#set setname = 'A'
#set nset    = 1

#set setname = 'A-G'
#set nset    = 7

set setname = 'A-U'
set nset    = 21

#*************************************************************************************

#set list_axis = "0 1 2 3 4 5 6" # 0(evsi) 1(abs-mmiss) 2(kfbchi) 3(kfbcl) 4(dzll) 5(cos-theta_B) 6(de)
#set list_axis = "2 3 4 6"
set list_axis = "6"
if( 1 == 1 ) then
  foreach axis ($list_axis)
    bsub -q ex ./exe_var.sh  $axis $stream $setname $nstream $nset
  end
endif

#set list_axis = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17" # 0-17
set list_axis = "2"
if( 1 == 0 ) then
  foreach axis ($list_axis)
    bsub -q ex ./exe_ksfw.sh  $axis $stream $setname $nstream $nset
  end
endif

