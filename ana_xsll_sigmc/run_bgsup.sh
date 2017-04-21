#!/bin/bash

#make;

for axis in 4 5 # 0 1 2 3 4 5
do
  bsub -q h ./exe_bgsup.sh $axis K-U # [fl_axis][setname]
done
