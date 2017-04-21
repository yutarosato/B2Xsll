#!/bin/csh -f

set list_set     = "A B C D E F G # H I J K L M N O P Q R S T U"
set list_xsid    = "1 2 3 4 5 6"
set outdir_merge = `dirname ${3}`/`basename ${3}`_merge/

mkdir -p $3;
mkdir -p ${outdir_merge};
     
foreach set($list_set)
  foreach xsid($list_xsid)
    (./bcs_lr_sig $1 $set $xsid 0 $2 $3 $4     >> log/log_bcs_lr_sig_xsid${xsid}_lep0_${4}.log ) >>& error.log
    (./bcs_lr_sig $1 $set $xsid 1 $2 $3 $4     >> log/log_bcs_lr_sig_xsid${xsid}_lep1_${4}.log ) >>& error.log
  end

  (./bcs_merge  $1 $set $3 ${outdir_merge}     >> log/log_bcs_merge_set${set}.log              ) >>& error.log
end

