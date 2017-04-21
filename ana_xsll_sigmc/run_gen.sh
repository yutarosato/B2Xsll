#!/bin/bash

make;

#setname="A-U" # A, AB, A-C
#nset="21"

setname="A-G" # A, AB, A-C
nset="7"

bsub -q ex ./exe_gen_M_xs.sh 2 1 $setname $nset # [fl_mode_ll][fl_klong][setname][nset]
bsub -q ex ./exe_gen_M_xs.sh 2 0 $setname $nset 
