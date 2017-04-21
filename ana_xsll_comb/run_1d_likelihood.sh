#!/bin/csh -f

#set stream  = '0'   # small dataset for debug
#set nstream = '1'   # small dataset for debug
#set setname = 'A'   # small dataset for debug
#set nset    = '1'   # small dataset for debug
set stream  = '0-2' # official dataset for significance
set nstream = '3'   # official dataset for significance
set setname = 'A-J' # official dataset for significance
set nset    = '10'  # official dataset for significance

##########################################################################################################

#<LR>          
#set tag     = '1_tot'  # {1,2,1ks,2ks}_tot
#set brname  = "lr${tag}"

#<NB>
set tag    = "madeksfw_vtxcl_fmiss1_tot"
set brname = "nb_${tag}"

bsub -q s ./exe_1d_likelihood.sh  $stream $setname $nstream $nset $tag $brname
##########################################################################################################

#<NB>(lep)
#set tag    = "orgksfw_vtxcl_fmiss1_tot"
#set brname = "nb_lep%d_${tag}"

#bsub -q s ./exe_1d_likelihood_lep.sh   $stream $setname $nstream $nset $tag $brname
#bsub -q s ./exe_1d_likelihood_lep5.sh  $stream $setname $nstream $nset $tag $brname
##########################################################################################################
