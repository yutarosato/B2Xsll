#!/bin/csh -f

#make;

#set setname = 'A-U'
set setname = 'A-G'
#set setname = 'A'

bsub -q ex ./exe_2d_q2_theta_eff.sh  1 $setname
bsub -q ex ./exe_2d_q2_theta_eff.sh  0 $setname

#bsub -q ex ./exe_2d_q2_theta_eff_transition10.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_transition10.sh  0 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_transition12.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_transition12.sh  0 $setname

#bsub -q ex ./exe_2d_q2_theta_eff_hadronization.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_hadronization.sh  0 $setname

#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kll_p.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kll_p.sh  0 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kll_m.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kll_m.sh  0 $setname

#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kstrll_p.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kstrll_p.sh  0 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kstrll_m.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_kstrll_m.sh  0 $setname

#bsub -q ex ./exe_2d_q2_theta_eff_fraction_xsll_p.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_xsll_p.sh  0 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_xsll_m.sh  1 $setname
#bsub -q ex ./exe_2d_q2_theta_eff_fraction_xsll_m.sh  0 $setname
