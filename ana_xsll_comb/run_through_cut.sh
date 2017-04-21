#!/bin/csh -f

#set  indir = "NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs/gMC_*s0[6-9]*.root"           #   gMC
#set outdir = "hbk_afb_calib/hbk_cut2/bkg_522/"

#set  indir = "NB_lep_calib/hbk_rd_cc/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs/RD_*.root"           #   RD
#set outdir = "hbk_afb_calib/hbk_cut3/rd_522/"

#set  indir = 'NB_lep_calib/hbk_c10flip/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_cut2/sig_c10flip_522/'

#set  indir = "../data/gmc2/hbk6/right/hbk_calib_cut3_nb/hbk_orgksfw_vtxcl_fmiss1/gMC*s0[0-5]*.root"           #   gMC2
#set outdir = "../data/gmc2/hbk6/right/hbk_calib_cut3_nb/hbk_orgksfw_vtxcl_fmiss1_lrnb/"

#set  indir = "hbk_bcs/sig/sig1_bb_522_merge/*.root"
#set outdir = "hbk_bcs/sig_lrnb/sig1_bb_522_merge/"

#set  indir = "../data/gmc_cc/hbk6/right/hbk_calib_cut5_nb/hbk_orgksfw_vtxcl_fmiss1_bb_bcs_522/*.root"
#set outdir = "hbk_afb_calib/hbk_cut5/cc_522/"

#set  indir = '../data/sigmc/hbk6/right/hbk_calib_cut2_large_nb/hbk_orgksfw_vtxcl_fmiss1_bb_bcs_merge_522/*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_cut2/sig_522_large/'

#set  indir = '../data/sigmc/hbk6/right/hbk_calib_cut2_c10flip_large_nb/hbk_orgksfw_vtxcl_fmiss1_bb_bcs_merge_522/*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_cut2/sig_c10flip_522_large/'

#set  indir = "NB_lep_calib/hbk_totq2/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs/gMC_*_s0[0-5]*.root"           #   gMC
#set outdir = "hbk_afb_calib/hbk_cut5/bkg_522/"

#set  indir = "NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/gMC*s0[0-5]*.root"           #   gMC
#set outdir = "tmp_gmc/"

#set  indir = 'NB_lep_calib/hbk_mb/hbk_orgksfw_vtxcl_fmiss1_mb465{,_c10flip}_bb_522_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_mb/sig_522_mb465/'

#set  indir = 'NB_lep_calib/hbk_fm/hbk_orgksfw_vtxcl_fmiss1_fm200{,_c10flip}_bb_522_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_fm/sig_522_fm200/'

#set  indir = 'NB_lep_calib/hbk_xsspin/hbk_orgksfw_vtxcl_fmiss1{,_c10flip}_bb_522_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_xsspin/sig_522/'

#set  indir = 'NB_lep_calib/hbk_lambdaone/hbk_orgksfw_vtxcl_fmiss1_lambdaone362{,_c10flip}_bb_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_lambdaone/sig_522_lambdaone362/'

#set  indir = 'tmp_hbk_forbcs_recmode/hbk_cut2/sig_522_lrnb_bb_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'tmp_hbk_forbcs_recmode/hbk_cut2/sig_522_lrnb_bb_bcs_merge_lrnb/'

#set  indir = "NB_lep_calib/hbk_cut8/hbk_orgksfw_vtxcl_fmiss1_bb_bcs/gMC_*_s0[0-5]*.root"           #   gMC
#set outdir = "hbk_afb_calib/hbk_cut8/bkg_522/"

#set  indir = "NB_lep_calib/hbk_rd_cut8/hbk_orgksfw_vtxcl_fmiss1_bb_bcs/RD_*.root"           #   RD
#set outdir = "hbk_afb_calib/hbk_cut8/rd_522/"

#set  indir = 'NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/sigMC*_set[A]*.root' # sigMC (bcs probability check)
#set outdir = 'hbk_for_bcsprob/hbk_cut2/sig_522/'

#set  indir = "NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/gMC_*_s0[0]*.root" # gMC (bcs probability check)
#set outdir = 'hbk_for_bcsprob/hbk_cut2/bkg_522/'

#set  indir = "NB_lep_calib/hbk_rd/hbk_orgksfw_vtxcl_fmiss1/RD_*.root" # RD (bcs probability check)
#set outdir = 'hbk_for_bcsprob/hbk_cut2/rd_522/'

#set  indir = 'NB_lep_calib/hbk_cut3/hbk_orgksfw_vtxcl_fmiss1_bb_522_bcs_merge/sigMC*_set[A-U]*.root' # sigMC
#set outdir = 'hbk_afb_calib/hbk_cut3/sig_522/'

set  indir = 'NB_lep_calib/hbk_ff6/hbk_orgksfw_vtxcl_fmiss1_ff6{,_c10flip}_bb_bcs_merge/sigMC*_set[A-N]*.root' # sigMC
set outdir = 'hbk_afb_calib/hbk_ff6/sig_522_ff6/'

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
    #bsub -q e ./exe_through_cut_lrnb.sh $f $outdir
  echo "./exe_through_cut_lrnb.sh $f $outdir" >> tmp.log
  end
else
    echo 'canceled'
endif
