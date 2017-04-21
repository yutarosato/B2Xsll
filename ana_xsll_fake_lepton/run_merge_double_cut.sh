#!/bin/csh -f

set  tag   = 'right' # {right,emu}
set indir  = "hbk_double_nb/hbk_orgksfw_vtxcl_fmiss1/"
set outdir = "hbk_double_nb/hbk_orgksfw_vtxcl_fmiss1_lrnb/"

#set list_exp = `awk '{print $1}' ~/script/nBB2.txt`
set list_exp = "13 15 17 21 31 35 39 51 65"
set list_stream = "0" # 0 1 2 3 4 5
set list_type = "rd"
#set list_type = "uds charm mixed charged"

mkdir -p $outdir
foreach exp ($list_exp)
   foreach stream ($list_stream)
      foreach type ($list_type)
         bsub -q s ./exe_merge_double_lrnb_cut2.sh $indir $outdir $exp $stream $type
      end
   end
end
