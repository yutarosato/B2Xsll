#!/bin/csh -f

set  indir = "~/ewp/ana/data/sigmc/hbk6/right/hbk_calib_cut2/*.root"
#set  indir = "~/ewp/ana/data/gmc/hbk6/right/hbk_calib_cut2/*.root"
set outdir = "tmp_hbk_cal_delh/"

set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute cal_delh(y/n) ?"
#set kakunin = $<
set  kakunin = 'y'
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
    bsub -q e ./exe_cal_delh.sh  $f $outdir
  end
else
    echo 'canceled'
endif
