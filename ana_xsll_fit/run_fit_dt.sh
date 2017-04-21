#!/bin/csh -f

#bsub -q s ./exe_Mbc_fit_dt.sh 15 ;


#bsub -q s ./exe_Mbc_fit_bin_dt_totFB.sh 15;
#exit;

#bsub -q s ./exe_Mbc_fit_bin_dt.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_mod.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod.sh 15 3;

bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 0;
bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 1;
bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 2;
bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_1_6q2.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_1_6q2.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_1_6q2.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_mod_PRL_1_6q2.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_kll.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_kll.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_kll.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_kll.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_kstrll.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_kstrll.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_kstrll.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_kstrll.sh 15 3;

#bsub -q s ./exe_Mbc_fit_bin_dt_xsll.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_xsll.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_xsll.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_xsll.sh 15 3;
exit

# [CF slope]
#foreach axis (`seq 0 1 15`)
#foreach axis (20 21)
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfslope.sh  15  0   $axis   # [sel_func][fl_q2][fl_axis]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfslope.sh  15  1   $axis   # [sel_func][fl_q2][fl_axis]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfslope.sh  15  2   $axis   # [sel_func][fl_q2][fl_axis]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfslope.sh  15  3   $axis   # [sel_func][fl_q2][fl_axis]
#end
#exit

# [signal shape, CF width, scf]
foreach cnt (`seq 1 1 100`)
#foreach cnt (`seq 101 1 200`)
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_sigshape.sh  15  0   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_sigshape.sh  15  1   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_sigshape.sh  15  2   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_sigshape.sh  15  3   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfwidth.sh  15  0   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfwidth.sh  15  1   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfwidth.sh  15  2   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfwidth.sh  15  3   $cnt   # [sel_func][fl_q2][fl_cnt]
  bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfwidth_beta.sh  15  1   $cnt   # [sel_func][fl_q2][fl_cnt]
  bsub -q s ./exe_Mbc_fit_bin_dt_syst_cfwidth_beta.sh  15  2   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf.sh  15  0   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf.sh  15  1   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf.sh  15  2   $cnt   # [sel_func][fl_q2][fl_cnt]
  #bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf.sh  15  3   $cnt   # [sel_func][fl_q2][fl_cnt]
end

# [scf_f]
#bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf_f.sh 15 0;
#bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf_f.sh 15 1;
#bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf_f.sh 15 2;
#bsub -q s ./exe_Mbc_fit_bin_dt_syst_scf_f.sh 15 3;


# [peaking B.G.]
#foreach q2 (`seq 0 1 3`)
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_charmonium_gmc.sh  15  $q2     # [sel_func][fl_q2]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_charmonium.sh      15  $q2  1  # [sel_func][fl_q2][fl_syst]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_charmonium.sh      15  $q2 -1  # [sel_func][fl_q2][fl_syst]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_double.sh          15  $q2  1  # [sel_func][fl_q2][fl_syst]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_double.sh          15  $q2 -1  # [sel_func][fl_q2][fl_syst]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_swap.sh            15  $q2  1  # [sel_func][fl_q2][fl_syst]
#  bsub -q s ./exe_Mbc_fit_bin_dt_syst_peak_swap.sh            15  $q2 -1  # [sel_func][fl_q2][fl_syst]
#end
