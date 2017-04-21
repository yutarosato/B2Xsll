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

#set stream  = '0'
#set nstream = '1'
#set setname = 'A-U'
#set nset    = '21'


############################[total and BIN]####################################
#foreach ratio (`seq 50 50 600`)
#foreach ratio (500)
   #bsub -q e ./exe_Mbc_fit_linearity.sh      $stream $setname $nstream $nset 15 $ratio; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#   bsub -q e ./exe_Mbc_fit_linearity_bin.sh  $stream $setname $nstream $nset 15 $ratio; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
   #bsub -q e ./exe_Mbc_fit_linearity_bin.sh  $stream $setname 1        $nset 15 $ratio; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#end

############################[BIN2(toyMC)]####################################
#foreach fl_cnt (`seq 1 1 9`)
#foreach ratio (`seq 10 5 300`)
#foreach fl_q2 (`seq 0 1 3`)
#foreach ratio (10)
#foreach fl_q2 (0)
#foreach fl_cnt (3)
#   bsub -q e ./exe_Mbc_fit_linearity_bin2_toy.sh  $stream $setname $nstream $nset 15 $ratio $fl_q2 $fl_cnt; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#end
#end
#end
  
############################[BIN2(ensemble test)]####################################
#[CAUTION] ROOT fit(iterative fit) should be skipped. 
foreach ratio (100)
#foreach fl_q2 (`seq 0 1 3`)
#foreach fl_q2 (2 3)
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  0 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  1 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  2 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  3 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  4 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  5 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]

  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  6 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  7 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  8 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
  #bsub -q e ./exe_Mbc_fit_linearity_bin2_ensemble.sh  9 $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#end
#end

############################[BIN3(toyMC)]####################################
#foreach fl_cnt (`seq 1 1 40`)
#foreach ratio (`seq -100 5 100`)
#foreach fl_q2 (`seq 0 1 3`)
#foreach ratio (-90)
#foreach fl_q2 (0)
#foreach fl_cnt (9)
#   bsub -q e ./exe_Mbc_fit_linearity_bin3_toy.sh  $stream $setname $nstream $nset 15 $ratio $fl_q2 $fl_cnt; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#end
#end
#end

############################[BIN3*obsolete]####################################
#foreach ratio (`seq -100 5 100`)
#foreach fl_q2 (`seq 0 1 3`)
#foreach ratio (-100 -5 -10 -20)
#foreach fl_q2 (1)
   #bsub -q e ./exe_Mbc_fit_linearity_bin3.sh  $stream $setname $nstream $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
   #bsub -q e ./exe_Mbc_fit_linearity_bin3.sh  $stream $setname 1        $nset 15 $ratio $fl_q2; # [stream][setname][nstream][nset] [sel_func] [incl_ratio]
#end
#end


#[SHAPE]
#foreach bin (`seq 0 1 8`)
#bsub -q e ./exe_Mbc_fit_shape.sh            8 $bin
#bsub -q e ./exe_Mbc_fit_shape_for_rd_emu.sh 8 $bin
#end

#foreach bin (`seq 0 1 14`)
#bsub -q e ./exe_Mbc_fit_shape.sh            14 $bin
#bsub -q e ./exe_Mbc_fit_shape_for_rd_emu.sh 14 $bin
#end
