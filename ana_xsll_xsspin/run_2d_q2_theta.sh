#!/bin/csh -f

#make;

set setname = 'A-U'
#set setname = 'A'

foreach xsid (0 1 2 3 4 5) # 0(all), 1(K), 2(K*), 3(Xs), 4(Xs-spin0), 5(Xs-spin1)
  bsub -q ex ./exe_2d_q2_theta_eff.sh  1  $xsid $setname
  bsub -q ex ./exe_2d_q2_theta_eff.sh  0  $xsid $setname
end


