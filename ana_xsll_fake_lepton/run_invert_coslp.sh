#!/bin/csh -f

set  indir = "../data/gmc_fl/hbk4/right/hbk_double_cut2/rd/9999/*.root"
set outdir = "../data/gmc_fl/hbk5/right/hbk_double_cut2/rd/9999/"

set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute calib(y/n) ?"
set kakunin = $<
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
    bsub -q s ./exe_invert_coslp.sh  $f $outdir
  end
else
    echo 'canceled'
endif
