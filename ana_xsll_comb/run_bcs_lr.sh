#!/bin/csh -f


#set tag    = '1_qq'  # {1,2,1ks,2ks}{_tot,_qq,_bb}
set tag    = $1  # {1,2,1ks,2ks}{_tot,_qq,_bb}
set brname = "lr${tag}"
#####################################################################

# <SIG>
set list_exp_sig = `awk '{print $1}' ~/script/nBB.txt`
set indir_sig    = "hbk_lr/sig/"
set outdir_sig   = "hbk_bcs/sig/sig${tag}/"

mkdir -p $outdir_sig
foreach exp($list_exp_sig)
    bsub -q s ./exe_bcs_lr_sig.sh            $exp $indir_sig $outdir_sig $brname
end

#####################################################################
# <BKG>
set list_exp_bkg = `awk '{print $1}' ~/script/nBB2.txt`
set list_stream  = "0 1 2 3 4 5" # 0-5
set indir_bkg    = "hbk_lr/bkg/"
set outdir_bkg   = "hbk_bcs/bkg/bkg${tag}/"

mkdir -p $outdir_bkg
foreach exp ($list_exp_bkg)
   foreach stream ($list_stream)
         bsub -q e ./exe_bcs_lr_bkg.sh           $exp $stream $indir_bkg $outdir_bkg $brname
   end
end


