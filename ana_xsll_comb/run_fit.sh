#!/bin/csh -f

#make;

#set stream  = '0-2' # evaluation
#set nstream = '3'   # evaluation
#set setname = 'A-J' # evaluation
#set nset    = '10'  # evaluation

set stream  = '0-5'
set nstream = '6'
set setname = 'A-U'
set nset    = '21'


bsub -q s ./exe_Mbc_fit.sh  $stream $setname $nstream $nset 15  # [stream][setname][nstream][nset] [sel_func]
#bsub -q s ./exe_Mbc_fit.sh  $stream $setname $nstream $nset 151

#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15 100 # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  90
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  80
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  70
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  60
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  50
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  40
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  30
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  20
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15  10
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15   0
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15 110
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15 120
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15 130
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15 140
#bsub -q s ./exe_Mbc_fit_linearity.sh  $stream $setname $nstream $nset 15 150
