#!/bin/csh -f

#make;

set stream  = '0-2' # evaluation
set nstream = '3'   # evaluation
set setname = 'A-J' # evaluation
set nset    = '10'  # evaluation

#set stream  = '0-5'
#set nstream = '6'
#set setname = 'A-J'
#set nset    = '10'


bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 1 0
bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 0 0

#bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 1 -1  # [stream][setname][nstream][nset] [fl_mode_ll][fl_q2]
#bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 0 -1
#bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 1  0
#bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 0  0
#bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 1  1
#bsub -q s ./exe_Mbc.sh  $stream $setname $nstream $nset 0  1

#./exe_Mbc_tmp.sh  $stream $setname $nstream $nset 0

if( 0 == 1 ) then
  foreach m (1 10    101 1001 110 1010    201 1101 210 1110    301 1201 310 1210    401 1301 410 1310)
    bsub -q s ./exe_Mbc_mode.sh  $stream $setname $nstream $nset 0  0 $m  #  [stream][setname][fl_modell][fl_q2][fl_mode_xs]
    bsub -q s ./exe_Mbc_mode.sh  $stream $setname $nstream $nset 1  0 $m
    bsub -q s ./exe_Mbc_mode.sh  $stream $setname $nstream $nset 0  1 $m
    bsub -q s ./exe_Mbc_mode.sh  $stream $setname $nstream $nset 1  1 $m
    bsub -q s ./exe_Mbc_mode.sh  $stream $setname $nstream $nset 0 -1 $m
    bsub -q s ./exe_Mbc_mode.sh  $stream $setname $nstream $nset 1 -1 $m
  end
endif

if( 0 == 1 ) then
   bsub -q s ./exe_Mbc_body.sh  $stream $setname $nstream $nset 1
   bsub -q s ./exe_Mbc_body.sh  $stream $setname $nstream $nset 2
   bsub -q s ./exe_Mbc_body.sh  $stream $setname $nstream $nset 3
   bsub -q s ./exe_Mbc_body.sh  $stream $setname $nstream $nset 4
   #bsub -q s ./exe_Mbc_body.sh  $stream $setname $nstream $nset 5
endif

