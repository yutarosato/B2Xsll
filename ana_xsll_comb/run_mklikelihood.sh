#!/bin/csh -f

set tag = "1" # 1, 2, 1ks, 2ks
# [SIG]
#set  indir = 'fmiss/hbk/sig/*_set[A-U]*'
#set outdir = "hbk_lr/sig/sig${tag}/"
# [BKG]
set  indir = 'fmiss/hbk/bkg/gMC_*_s0[0-5]*'
set outdir = "hbk_lr/bkg/bkg${tag}/"

set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute mklikelihood(y/n) ?"
set kakunin = $<
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
     bsub -q s ./exe_mklikelihood.sh $f $outdir $tag
  end
else
    echo 'canceled'
endif

#./exe_mklikelihood.sh "plot" 1 1
