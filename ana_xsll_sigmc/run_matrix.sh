#!/bin/bash

#make;

################################################################################
# count the number of generated events
################################################################################
if [ 1 == 0 ]; then
  for set in A B C D E F G H I J K L M N O P Q R S T U
    do
     #[all]
       bsub -q e ./exe_entry.sh 0 0 $set # [fl_mode_ll] [fl_xsid] [setname]
       bsub -q e ./exe_entry.sh 1 0 $set
     #[K]
       bsub -q e ./exe_entry.sh 0 1 $set
       bsub -q e ./exe_entry.sh 1 1 $set
     #[K*]
       bsub -q e ./exe_entry.sh 0 2 $set
       bsub -q e ./exe_entry.sh 1 2 $set
     #[Xs]
       bsub -q e ./exe_entry.sh 0 3 $set
       bsub -q e ./exe_entry.sh 1 3 $set
  done
fi


################################################################################
#  count the number of reconstructed events
################################################################################       
if [ 1 == 1 ]; then
  for set in  A B C D E F G H I J K L M N O P Q R S T U 
    do
      for xs in 1 10    101 1001 110 1010    201 1101 210 1110    301 1201 310 1210    401 1301 410 1310
        do
          #[all]
	    #bsub -q e ./exe_all.sh      $xs 0 0 $set # [fl_mode_xs] [fl_mode_ll] [fl_xsid] [setname]
	    #bsub -q e ./exe_all.sh      $xs 1 0 $set
          #[K]
	    #bsub -q e ./exe_all.sh      $xs 0 1 $set
	    #bsub -q e ./exe_all.sh      $xs 1 1 $set
          #[K*]
	    #bsub -q e ./exe_all.sh      $xs 0 2 $set
	    #bsub -q e ./exe_all.sh      $xs 1 2 $set
          #[Xs]
   	    bsub -q s ./exe_all.sh      $xs 0 3 $set
	    bsub -q s ./exe_all.sh      $xs 1 3 $set
       done
  done
fi
