#!/bin/csh -f

#set  indir = "~/ewp/ana/data/sigmc/hbk6/right/hbk_org_c10flip_large/sigMC_*_set[A-U]*.root"
#set outdir = "~/ewp/ana/data/sigmc/hbk6/right/hbk_calib_c10flip_large/"
#set  indir = "~/ewp/ana/data/gmc/hbk6/right/hbk_org/mixed/9999/gMC_*_*s0[6-9]*.root"
#set outdir = "~/ewp/ana/data/gmc/hbk6/right/hbk_calib/mixed/9999/"
#set  indir = "~/ewp/ana/data/gmc/hbk6/right/hbk_org/rd/9999/RD_*.root"
#set outdir = "~/ewp/ana/data/gmc/hbk6/right/hbk_calib/rd/9999/"
#set  indir = "~/ewp/ana/data/sigmc/hbk6/right/hbk_org_c10flip/sigMC_*_set[A-U]*.root"
#set outdir = "~/ewp/ana/data/sigmc/hbk6/right/hbk_calib_c10flip/"
#set  indir = "~/ewp/ana/data/gmc_cc/hbk6/right/hbk_org/CC_mixedjpsi_*.root"
#set outdir = "~/ewp/ana/data/gmc_cc/hbk6/right/hbk_calib/"

#set  indir = "~/ewp/ana/data/gmc/hbk6/emu/hbk_org/${1}/9999/gMC_*s0[5]*.root"
#set outdir = "~/ewp/ana/data/gmc/hbk6/emu/hbk_calib/${1}/9999/"
#set  indir = "~/ewp/ana/data/gmc/hbk6/emu/hbk_org/rd/9999/RD_*.root"
#set outdir = "~/ewp/ana/data/gmc/hbk6/emu/hbk_calib/rd/9999/"

#set  indir = "~/ewp/ana/data/sigmc_xsspin/hbk6/right/hbk_org_c10flip/*.root"
#set outdir = "~/ewp/ana/data/sigmc_xsspin/hbk6/right/hbk_calib_c10flip/"

#set  indir = "~/ewp/ana/data/sigmc_lambdaone/hbk6/right/hbk_org_lambdaone429/*.root"
#set outdir = "~/ewp/ana/data/sigmc_lambdaone/hbk6/right/hbk_calib_lambdaone429/"

#set  indir = "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd_cc/hbk_orgksfw_vtxcl_fmiss1/*.root"
#set outdir = "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd_cc/hbk_orgksfw_vtxcl_fmiss1_through_calib_code/"

set  indir = "~/ewp/ana/data/sigmc_ff/hbk6/right/hbk_org_ff6_c10flip/*set[A-G].root"
set outdir = "~/ewp/ana/data/sigmc_ff/hbk6/right/hbk_calib_ff6/"

set count=0

foreach f ($indir)
     @ count++
end
 echo $count "files exists"

echo "execute calib(y/n) ?"
#set kakunin = $<
set  kakunin = 'y'
@ count=0
if( $kakunin == 'y') then
  mkdir -p $outdir
  foreach f ($indir)
    @ count++
    #bsub -q e ./exe_calib.sh  $f $outdir
    echo "./exe_calib.sh  $f $outdir" >> tmp.log
  end
else
    echo 'canceled'
endif
