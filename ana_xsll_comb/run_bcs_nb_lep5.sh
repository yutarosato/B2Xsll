#!/bin/csh -f


set tag    = $1
set type   = $2
#set tag    = "orgksfw_vtxcl_fmiss1"
#set type   = "bb"      # tot, qq, bb
set brname = "nb_lep%d_${tag}_${type}"
#####################################################################
set indir    = "NB_lep5/hbk/hbk_${tag}/"
set outdir   = "NB_lep5/hbk/hbk_${tag}_${type}_bcs/"
#set outdir   = "NB_lep5/hbk_5body/hbk_${tag}_${type}_bcs/" # remove K4pi modes
#####################################################################

# <SIG>
set list_exp_sig = `awk '{print $1}' ~/script/nBB.txt`

foreach exp($list_exp_sig)
    bsub -q s ./exe_bcs_lr_lep_sig.sh            $exp $indir $outdir $brname
end

#####################################################################
# <BKG>
set list_exp_bkg = `awk '{print $1}' ~/script/nBB2.txt`
set list_stream  = "0 1 2 3 4 5" # 0-5

foreach exp ($list_exp_bkg)
   foreach stream ($list_stream)
     bsub -q e ./exe_bcs_lr_lep_bkg.sh           $exp $stream $indir $outdir $brname
   end
end
