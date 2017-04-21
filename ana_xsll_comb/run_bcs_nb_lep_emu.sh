#!/bin/csh -f

set tag    = $1
set type   = $2
#set tag    = "orgksfw_vtxcl_fmiss1"
#set type   = "tot"      # tot, qq, bb
set brname = "nb_lep%d_${tag}_${type}"
#####################################################################
set list_exp_bkg = `awk '{print $1}' ~/script/nBB2.txt`
# <BKG>
set list_stream  = "0 1 2 3 4 5" # 0-5
set indir_bkg    = "NB_lep_calib/hbk_emu/hbk_${tag}/"
set outdir_bkg   = "NB_lep_calib/hbk_emu/hbk_${tag}_${type}_bcs/"

#foreach exp ($list_exp_bkg)
#   foreach stream ($list_stream)
#         bsub -q e ./exe_bcs_lr_lep_bkg_emu.sh           $exp $stream $indir_bkg $outdir_bkg $brname
#   end
#end

# <RD>
set indir_bkg    = "NB_lep_calib/hbk_rd_emu/hbk_${tag}/"
set outdir_bkg   = "NB_lep_calib/hbk_rd_emu/hbk_${tag}_${type}_bcs/"

foreach exp ($list_exp_bkg)
    bsub -q e ./exe_bcs_lr_lep_rd_emu.sh           $exp $indir_bkg $outdir_bkg $brname
end
