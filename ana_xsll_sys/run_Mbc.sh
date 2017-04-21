#!/bin/csh -f

#set setname = 'A'
#set setname = 'A-J'
set setname = 'A-U'


#foreach xs (1 10    101 1001 110 1010    201 1101 210 1110    301 1201 310 1210    401 1301 410 1310)
#   bsub -q e ./exe_Mbc_mode.sh $xs 0 $setname
#   bsub -q e ./exe_Mbc_mode.sh $xs 1 $setname
#end



bsub -q e ./exe_Mbc_Ks_pi0.sh -1000 0 $setname # [fl_mode][fl_mc][setname]
bsub -q e ./exe_Mbc_Ks_pi0.sh -1000 1 $setname
bsub -q e ./exe_Mbc_Ks_pi0.sh  1000 0 $setname
bsub -q e ./exe_Mbc_Ks_pi0.sh  1000 1 $setname
bsub -q e ./exe_Mbc_Ks_pi0.sh   -10 0 $setname
bsub -q e ./exe_Mbc_Ks_pi0.sh   -10 1 $setname
bsub -q e ./exe_Mbc_Ks_pi0.sh    10 0 $setname
bsub -q e ./exe_Mbc_Ks_pi0.sh    10 1 $setname
