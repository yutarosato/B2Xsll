#!/bin/csh -f

#set  indir = "NB_lep_calib/hbk_emu/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs/gMC*s0[0-5]*.root"           #   gMC(emu)
#set outdir = "hbk_afb_calib/hbk_cut2/bkg_emu_522/"
set  indir = "NB_lep_calib/hbk_rd_emu/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs/RD_*.root"           #   RD(emu)
set outdir = "hbk_afb_calib/hbk_cut2/rd_emu_522/"

set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute through_cut(y/n) ?"
set kakunin = $<
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
    bsub -q e ./exe_through_cut_lrnb_emu.sh $f $outdir 0
    bsub -q e ./exe_through_cut_lrnb_emu.sh $f $outdir 1
  end
else
    echo 'canceled'
endif
