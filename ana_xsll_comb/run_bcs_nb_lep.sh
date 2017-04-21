#!/bin/csh -f


set tag    = $1
set type   = $2
#set tag    = "orgksfw_vtxcl_fmiss1"
#set type   = "bb"      # tot, qq, bb
set brname = "nb_lep%d_${tag}_${type}"
#####################################################################
#set indir    = "../ana_xsll_fake_lepton/hbk_single_jpsi/"
#set outdir   = "../ana_xsll_fake_lepton/hbk_single_jpsi_bcs/"
#set indir    = "NB_lep_calib_test/hbk/hbk_${tag}/"
#set outdir   = "NB_lep_calib_test/hbk/hbk_${tag}_${type}_bcs/"
#set indir    = "NB_lep_calib/hbk_mb/hbk_${tag}_mb495/"
#set outdir   = "NB_lep_calib/hbk_mb/hbk_${tag}_mb495_${type}_bcs/"
#set indir    = "NB_lep_calib/hbk_xsspin/hbk_${tag}/"
#set outdir   = "NB_lep_calib/hbk_xsspin/hbk_${tag}_${type}_bcs/"
#set indir    = "NB_lep_calib/hbk_lambdaone/hbk_${tag}_lambdaone362_c10flip/"
#set outdir   = "NB_lep_calib/hbk_lambdaone/hbk_${tag}_lambdaone362_c10flip_${type}_bcs/"
#set indir    = "NB_lep_calib/hbk_rd_cut8/hbk_${tag}/"
#set outdir   = "NB_lep_calib/hbk_rd_cut8/hbk_${tag}_${type}_bcs/"

#set indir    = "hbk_for_bcsprob/hbk_cut2/sig_522/"
#set outdir   = "hbk_for_bcsprob/hbk_cut2/sig_522_${type}_bcs/"
#set indir    = "hbk_for_bcsprob/hbk_cut2/bkg_522/"
#set outdir   = "hbk_for_bcsprob/hbk_cut2/bkg_522_${type}_bcs/"
#set indir    = "hbk_for_bcsprob/hbk_cut2/rd_522/"
#set outdir   = "hbk_for_bcsprob/hbk_cut2/rd_522_${type}_bcs/"

#set indir    = "NB_lep_calib/hbk_cut3_c10flip/hbk_${tag}/"
#set outdir   = "NB_lep_calib/hbk_cut3_c10flip/hbk_${tag}_${type}_bcs/"

set indir    = "NB_lep_calib/hbk_ff6/hbk_${tag}_ff6/"
set outdir   = "NB_lep_calib/hbk_ff6/hbk_${tag}_ff6_${type}_bcs/"
#set indir    = "NB_lep_calib/hbk_ff6/hbk_${tag}_ff6_c10flip/"
#set outdir   = "NB_lep_calib/hbk_ff6/hbk_${tag}_ff6_c10flip_${type}_bcs/"
#####################################################################

# <SIG>
set list_exp_sig = `awk '{print $1}' ~/script/nBB.txt`

foreach exp($list_exp_sig)
  bsub -q e ./exe_bcs_lr_lep_sig.sh            $exp $indir $outdir $brname
end
exit
#####################################################################
# <BKG>
set list_exp_bkg = `awk '{print $1}' ~/script/nBB2.txt`
set list_stream  = "0" # 0-5
#set list_stream  = "0 1 2 3 4 5" # 0-5
#set list_stream  = "6 7 8 9" # 0-5

#foreach exp ($list_exp_bkg)
#   foreach stream ($list_stream)
#     bsub -q e ./exe_bcs_lr_lep_bkg.sh           $exp $stream $indir $outdir $brname
#   end
#end

#####################################################################
# <RD>
foreach exp ($list_exp_bkg)
     bsub -q e ./exe_bcs_lr_lep_rd.sh           $exp $indir $outdir $brname
end
exit
#####################################################################
# <CC MC>

#set indir    = "../data/gmc_cc/hbk5/right/hbk_calib_cut2_nb/hbk_${tag}/*"
#set outdir   = "../data/gmc_cc/hbk5/right/hbk_calib_cut2_nb/hbk_${tag}_${type}_bcs/"

set indir    = "../data/gmc_cc/hbk6/right/hbk_calib_cut5_nb/hbk_${tag}/*"
set outdir   = "../data/gmc_cc/hbk6/right/hbk_calib_cut5_nb/hbk_${tag}_${type}_bcs/"

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
    bsub -q e ./exe_bcs_lr_lep_bkg_cc.sh  $f $outdir $brname
  end
else
    echo 'canceled'
endif
