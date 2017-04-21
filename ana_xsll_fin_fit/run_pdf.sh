#!/bin/csh -f


# fl_sb should be checked in draws_1d_Mbc_peak_{double,swap}?.cpp

make;

bsub -q ex ./exe_Mbc_peak_cc2.sh     0-5 6
exit

bsub -q ex ./exe_Mbc_self_cf2.sh     A-U
bsub -q ex ./exe_Mbc_peak_cc2.sh     0 100
bsub -q ex ./exe_Mbc_peak_double2.sh 0   1
bsub -q ex ./exe_Mbc_peak_swap2.sh   0   1
