#!/bin/csh -f

set  tag   = 'right' # {right,emu} # if emu is selected, outfile-name should be changed in merge_cut2.sh
#set  indir = "~/ewp/ana/data/gmc/hbk6/${tag}/hbk_org/"
#set  indir = "~/ewp/ana/data/gmc/hbk6/${tag}/hbk_calib/"
set  indir = "~/ewp/ana/data/gmc/hbk6/${tag}/hbk_org/"

set list_exp = `awk '{print $1}' ~/script/nBB2.txt`
#set list_stream = "6 7 8 9" # 0 1 2 3 4 5
set list_stream = "0 1 2 3 4 5" # 0 1 2 3 4 5
#set list_stream = "0" # 0 1 2 3 4 5
#set list_type = "mixed charged" # uds charm mixed charged rd
set list_type = "uds charm mixed charged" # uds charm mixed charged rd

#foreach exp ($list_exp)
#   foreach stream ($list_stream)
#      foreach type ($list_type)
         #bsub -q e ./exe_merge_cut.sh  $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut/"  $exp $stream $type
         #bsub -q e ./exe_merge_cut2.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut2/" $exp $stream $type
         #bsub -q e ./exe_merge_cut3.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut3/" $exp $stream $type
         #bsub -q e ./exe_merge_cut4.sh $indir "${HOME}/ewp/ana/data/gmc2/hbk6/${tag}/hbk_cut4/" $exp $stream $type
         #bsub -q e ./exe_merge_cut2.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_calib_cut2/" $exp $stream $type
         #bsub -q e ./exe_merge_cut3.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_calib_cut3/" $exp $stream $type
         #bsub -q e ./exe_merge_cut5.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_calib_cut5/" $exp $stream $type
         #bsub -q e ./exe_merge_cut7.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_calib_cut7/" $exp $stream $type
         #bsub -q e ./exe_merge_cut8.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_calib_cut8/" $exp $stream $type
#      end
#   end
#end


foreach exp ($list_exp)
    #bsub -q e ./exe_merge_cut.sh  $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut/"  $exp 999 rd #  [exp] [stream]
    #bsub -q e ./exe_merge_cut2.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut2/" $exp 999 rd
    #bsub -q e ./exe_merge_cut3.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut3/" $exp 999 rd
    bsub -q e ./exe_merge_cut8.sh $indir "${HOME}/ewp/ana/data/gmc/hbk6/${tag}/hbk_cut8/" $exp 999 rd
    #bsub -q e ./exe_merge_cut4.sh $indir "${HOME}/ewp/ana/data/gmc2/hbk6/${tag}/hbk_cut4/" $exp 999 rd
end
