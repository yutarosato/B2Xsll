#!/bin/csh -f

set TYPE   = $1
#set TYPE   = "uds" # uds ,charm, mixed, charged
set indir  = "~/ewp/ana/data/gmc_fl/hbk4/right/hbk_double_org/${TYPE}/9999/*.root"
set outdir = "~/ewp/ana/data/gmc_fl/hbk4/right/hbk_double_cut2/${TYPE}/9999/"



set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute through_cut(y/n) ?"
#set kakunin = $<
set kakunin = 'y'
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
  bsub -q e ./exe_through_double_cut2.sh $f $outdir
  end
else
    echo 'canceled'
endif
