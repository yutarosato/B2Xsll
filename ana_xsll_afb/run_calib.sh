#!/bin/csh -f

#make;

#set stream  = '0'
#set nstream = '1'
     
#set stream  = '3-5' # pdf
#set nstream = '3'   # pdf
     
#set stream  = '0-2' # evaluation
#set nstream = '3'   # evaluation
     
set stream  = '0-5' # all data
set nstream = '6'   # all data

#*************************************************************************************

#bsub -q s ./exe_dzll3d.sh       $stream $nstream  #          [stream] [nstream]
#bsub -q s ./exe_deltaE.sh     1 $stream $nstream  # [fl_pi0] [stream] [nstream]
#bsub -q s ./exe_deltaE.sh     0 $stream $nstream  # [fl_pi0] [stream] [nstream]

#*************************************************************************************

#bsub -q s ./exe_dzll3d.sh       0-5 6
#bsub -q s ./exe_bvtxchi2.sh     0-5 6

#bsub -q s ./exe_deltaE.sh     1 0-5 6
#bsub -q s ./exe_deltaE.sh     0 0-5 6

bsub -q s ./exe_bvtxchi2_scan.sh      0 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh      5 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     10 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     15 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     20 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     25 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     30 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     35 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     40 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     45 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     50 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     55 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     60 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     65 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     70 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     75 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     80 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     85 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     90 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh     95 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    100 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    105 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    110 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    115 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    120 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    125 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    130 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    135 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    140 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    145 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    150 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    155 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    160 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    165 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    170 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    175 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    180 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    185 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    190 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    195 0-5 6
bsub -q s ./exe_bvtxchi2_scan.sh    200 0-5 6
