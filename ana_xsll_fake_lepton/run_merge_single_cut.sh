#!/bin/csh -f

set  tag   = 'right' # {right,emu}
set indir  = "hbk_single_nb/hbk_orgksfw_vtxcl_fmiss1/"
set outdir = "hbk_single_nb/hbk_orgksfw_vtxcl_fmiss1_lrnb/"

set list_exp = `awk '{print $1}' ~/script/nBB2.txt`
set list_stream = "0" # 0 1 2 3 4 5
set list_type = "uds charm mixed charged"

mkdir -p $outdir
foreach exp ($list_exp)
   foreach stream ($list_stream)
      foreach type ($list_type)
         bsub -q s ./exe_merge_single_lrnb_cut2.sh $indir $outdir $exp $stream $type
      end
   end
end
