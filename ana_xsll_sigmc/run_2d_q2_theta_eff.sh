#!/bin/csh -f

#make;

#set setname = 'A'
#set setname = 'A-C'
#set setname = 'A-G'
#set setname = 'A-N'
set setname = 'A-U'

foreach xsid (0 1 2 3) # 0(all), 1(K), 2(K*), 3(Xs)
  #bsub -q e ./exe_2d_q2_theta_eff.sh       1  $xsid $setname # [fl_mode_ll] [setname]
  #bsub -q e ./exe_2d_q2_theta_eff.sh       0  $xsid $setname
  #bsub -q e ./exe_2d_q2_theta_eff_cut2.sh  1  $xsid $setname
  #bsub -q e ./exe_2d_q2_theta_eff_cut2.sh  0  $xsid $setname
  bsub -q e ./exe_2d_q2_theta_eff_cut4.sh  1  $xsid $setname
  bsub -q e ./exe_2d_q2_theta_eff_cut4.sh  0  $xsid $setname
end


