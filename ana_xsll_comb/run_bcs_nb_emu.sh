#!/bin/csh -f

#set tag    = $1
#set type    = $1
set tag    = "orgksfw_vtxcl_fmiss1"
set type   = "tot"      # tot, qq, bb
set brname = "nb_${tag}_${type}"
#####################################################################
# <BKG>
set list_exp_bkg = `awk '{print $1}' ~/script/nBB2.txt`
set list_stream  = "0 1 2 3 4 5" # 0-5
set indir_bkg    = "NB/hbk/bkg_emu/"
set outdir_bkg   = "NB/hbk/bkg_emu_${tag}_${type}_bcs/"

foreach exp ($list_exp_bkg)
  foreach stream ($list_stream)
     bsub -q e ./exe_bcs_lr_bkg_emu.sh           $exp $stream $indir_bkg $outdir_bkg $brname
  end
end
