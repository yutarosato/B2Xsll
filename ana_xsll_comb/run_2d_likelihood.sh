#!/bin/csh -f

#set stream  = '0'   # small dataset for debug
#set nstream = '1'   # small dataset for debug
#set setname = 'A'   # small dataset for debug
#set nset    = '1'   # small dataset for debug
set stream  = '0-2' # official dataset for significance
set nstream = '3'   # official dataset for significance
set setname = 'A-J' # official dataset for significance
set nset    = '10'  # official dataset for significance
#set stream  = '0-5' # official dataset for ALL
#set nstream = '6'   # official dataset for ALL
#set setname = 'A-U' # official dataset for ALL
#set nset    = '21'  # official dataset for ALL
##########################################################################################################
     
# LR     
set tag     = '1'  # {1,2,1ks,2ks}
set brname  = "lr${tag}"

# NB
#set tag    = "orgksfw_vtxcl_fmiss1"
#set brname = "nb_${tag}"

#bsub -q e ./exe_2d_likelihood.sh  $stream $setname $nstream $nset $tag "qq" $brname
#bsub -q e ./exe_2d_likelihood.sh  $stream $setname $nstream $nset $tag "bb" $brname
##########################################################################################################

# NB(lep)
set tag    = "orgksfw_vtxcl_fmiss1"
set brname = "nb_lep%d_${tag}"

#bsub -q e ./exe_2d_likelihood_lep.sh   $stream $setname $nstream $nset $tag "qq" $brname
bsub -q e ./exe_2d_likelihood_lep.sh   $stream $setname $nstream $nset $tag "bb" $brname
#bsub -q e ./exe_2d_likelihood_lep5.sh  $stream $setname $nstream $nset $tag "qq" $brname
#bsub -q e ./exe_2d_likelihood_lep5.sh  $stream $setname $nstream $nset $tag "bb" $brname
##########################################################################################################
