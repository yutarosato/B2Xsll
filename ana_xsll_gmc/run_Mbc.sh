#!/bin/csh -f

#make;

#set stream  = '0'   # pdf
#set nstream = '1'   # pdf
    
#set stream  = '3-5' # pdf
#set nstream = '3'   # pdf
     
#set stream  = '0-2' # evaluation
#set nstream = '3'   # evaluation
     
set stream  = '0-5' # all data
set nstream = '6'   # all data

#bsub -q s ./exe_Mbc_type.sh            $stream $nstream

#bsub -q s ./exe_Mbc_bb_lepself1.sh     $stream $nstream
#bsub -q s ./exe_Mbc_bb_lepself2.sh     $stream $nstream

#bsub -q s ./exe_Mbc_bb_smid1.sh      $stream $nstream
bsub -q s ./exe_Mbc_bb_smid2.sh      $stream $nstream

#bsub -q s ./exe_Mbc_bb_dmid1.sh    $stream $nstream
#bsub -q s ./exe_Mbc_bb_dmid2.sh    $stream $nstream

#bsub -q s ./exe_Mbc_bb_peak1.sh      $stream $nstream
#bsub -q s ./exe_Mbc_bb_peak2.sh      $stream $nstream

#bsub -q s ./exe_Mbc_bb_peak_tf1.sh      $stream $nstream
#bsub -q s ./exe_Mbc_bb_peak_tf2.sh      $stream $nstream

#bsub -q s ./exe_Mbc_bb_rest_peak2.sh      $stream $nstream
