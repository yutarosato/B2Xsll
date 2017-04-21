#!/bin/csh -f

#set indir = "../data/sigmc/hbk6/right/hbk_calib/*_set[A-U]*.root"
#set indir = "../data/sigmc/hbk6/right/hbk_calib_c10flip/*_set[A-U]*.root"
#set indir = "../data/sigmc/hbk6/right/hbk_calib_c10flip_large/*_set[A-U]*.root"
#set indir = "../data/sigmc/hbk6/right/hbk_calib_c10flip_large/*.root"
#set indir = "../data/gmc_cc/hbk6/right/hbk_calib/CC_mixedjpsi_*.root"
#set indir = "../data/sigmc/hbk6/right/hbk_org_c10sym_all/*_set?.root"

#set indir = "../data/sigmc_fm/hbk6/right/hbk_calib/*_set[A-U]*.root"
#set indir = "../data/sigmc_xsjpsi/hbk6/right/hbk_calib/*_set[A-J]*.root"
#set indir = "../data/sigmc_mb/hbk6/right/hbk_calib/*_set[A-U]*.root"

#set indir = "../data/sigmc_xsspin/hbk6/right/hbk_calib_c10flip/*_set[A-U]*.root"

#set indir = "../data/sigmc_lambdaone/hbk6/right/hbk_calib_lambdaone429/*_set[A-U]*.root"
set indir = "../data/sigmc_ff/hbk6/right/hbk_calib_ff6/*_set[A-N]*.root"

set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute through_cut(y/n) ?"
set kakunin = $<
@ count=0
if( $kakunin == 'y') then
  foreach f ($indir)
    @ count++
  #bsub -q e ./exe_through_cut.sh  $f '../data/sigmc/hbk6/right/hbk_cut/'
  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc/hbk6/right/hbk_cut2/'
  #bsub -q e ./exe_through_cut3.sh $f '../data/sigmc/hbk6/right/hbk_cut3/'
  #bsub -q e ./exe_through_cut4.sh $f '../data/sigmc/hbk6/right/hbk_cut4/'

  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut2_c10flip_large/'
  #bsub -q e ./exe_through_cut4.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut4_c10flip/'

  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc/hbk6/right/hbk_cut2/'
  #bsub -q e ./exe_through_cut3.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut3/'

  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc/hbk6/right/hbk_org_c10sym_all_tmp/'
  #bsub -q e ./exe_through_cut2.sh $f '../data/gmc_cc/hbk6/right/hbk_calib_cut2/'
  #bsub -q e ./exe_through_cut3.sh $f '../data/gmc_cc/hbk6/right/hbk_calib_cut3/'
  #bsub -q e ./exe_through_cut5.sh $f '../data/gmc_cc/hbk6/right/hbk_calib_cut5/'

  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc_fm/hbk6/right/hbk_calib_cut2/'
  #bsub -q e ./exe_through_cut3.sh $f '../data/sigmc_xsjpsi/hbk6/right/hbk_calib_cut3/'
  #bsub -q e ./exe_through_cut.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut/'

  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc_lambdaone/hbk6/right/hbk_calib_lambdaone429_cut2/'

  #bsub -q e ./exe_through_cut6.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut6/'
  #bsub -q e ./exe_through_cut6.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut6_c10flip/'

  #bsub -q e ./exe_through_cut2.sh $f '../data/sigmc_xsspin/hbk6/right/hbk_calib_cut2_c10flip/'
  #bsub -q e ./exe_through_cut7.sh $f '../data/sigmc/hbk6/right/hbk_calib_cut7/'
  #bsub -q e ./exe_through_cut2.sh $f 'tmp/'

  #echo "./exe_through_cut2.sh $f ../data/sigmc/hbk6/right/hbk_calib_cut2/" >> tmp.log
  echo "./exe_through_cut2.sh $f ../data/sigmc_ff/hbk6/right/hbk_calib_ff6_cut2/" >> tmp.log
  end
else
    echo 'canceled'
endif
