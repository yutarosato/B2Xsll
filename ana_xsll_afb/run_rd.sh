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

set setname = 'A-U'
set nset    = 21

#set setname = 'A'
#set nset    = 1

#*************************************************************************************

#bsub -q e ./exe_2d_q2_theta_gmc1.sh     $stream $nstream 1  # [stream] [nstream] [fl_mode_ll]
#bsub -q e ./exe_2d_q2_theta_gmc1.sh     $stream $nstream 0
#bsub -q e ./exe_2d_q2_theta_gmc2.sh     $stream $nstream 1
#bsub -q e ./exe_2d_q2_theta_gmc2.sh     $stream $nstream 0

#bsub -q e ./exe_2d_q2_theta_rd1.sh  $stream $nstream 1
#bsub -q e ./exe_2d_q2_theta_rd1.sh  $stream $nstream 0
#bsub -q e ./exe_2d_q2_theta_rd2.sh  $stream $nstream 1
#bsub -q e ./exe_2d_q2_theta_rd2.sh  $stream $nstream 0
     
#bsub -q e ./exe_2d_nb_lep_rd.sh  $stream $nstream 1  # [stream] [nstream] [fl_mode_ll]
#bsub -q e ./exe_2d_nb_lep_rd.sh  $stream $nstream 0

#bsub -q e ./exe_Mbc.sh       $stream $setname $nstream $nset
bsub -q sx ./exe_Mbc_sep_peak_scf.sh       $stream $setname $nstream $nset
#bsub -q e ./exe_Mbc_jpsi.sh  $stream $nstream
#bsub -q e ./exe_Mbc_emu.sh   $stream $nstream 0 15
#bsub -q e ./exe_Mbc_emu.sh   $stream $nstream 0 151
#bsub -q e ./exe_Mbc_emu.sh   $stream $nstream 0 50
#bsub -q e ./exe_Mbc_emu.sh   $stream $nstream 0 51

#bsub -q e ./exe_Mbc_body.sh   $stream $setname $nstream $nset 1
#bsub -q e ./exe_Mbc_body.sh   $stream $setname $nstream $nset 2
#bsub -q e ./exe_Mbc_body.sh   $stream $setname $nstream $nset 3
#bsub -q e ./exe_Mbc_body.sh   $stream $setname $nstream $nset 4
#bsub -q e ./exe_Mbc_body.sh   $stream $setname $nstream $nset 5

#set list_axis = "0 1 2 3 4 5 6" # 0(evsi) 1(abs-mmiss) 2(kfbchi) 3(kfbcl) 4(dzll) 5(cos-theta_B) 6(de)
set list_axis = "2 3 4 6"
if( 1 == 0 ) then
  foreach axis ($list_axis)
    bsub -q e ./exe_var.sh  $axis $stream $nstream
  end
endif

set list_axis = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17" # 0-17
if( 0 == 1 ) then
  foreach axis ($list_axis)
    bsub -q e ./exe_ksfw.sh  $axis $stream $nstream
  end
endif

