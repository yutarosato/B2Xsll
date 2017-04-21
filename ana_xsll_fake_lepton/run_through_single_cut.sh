#!/bin/csh -f

set TYPE   = $1
#set TYPE   = "uds" # uds ,charm, mixed, charged
set indir  = "~/ewp/ana/data/gmc_fl/hbk5/right/hbk_single_org/${TYPE}/9999/*.root"
set outdir = "~/ewp/ana/data/gmc_fl/hbk5/right/hbk_single_cut2/${TYPE}/9999/"



set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute through_cut(y/n) ?"
set kakunin = 'y'
#set kakunin = $<
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
    bsub -q s ./exe_through_single_cut2.sh $f $outdir
  end
else
    echo 'canceled'
endif
