#!/bin/csh -f

#set list_set     = "A B C D E F G"
#set list_set     = "H I J K L M N"
#set list_set     = "A B C D E F G H I J"
set list_set     = "A B C D E F G H I J K L M N O P Q R S T U"
#set list_set     = "AA       BB       CC       DD       EE       FF       GG       HH       II       JJ       KK       LL       MM       NN       OO       PP       QQ       RR       SS       TT       UU"
#set list_set     = "AAA      BBB      CCC      DDD      EEE      FFF      GGG      HHH      III      JJJ      KKK      LLL      MMM      NNN      OOO      PPP      QQQ      RRR      SSS      TTT      UUU"
#set list_set     = "AAAA     BBBB     CCCC     DDDD     EEEE     FFFF     GGGG     HHHH     IIII     JJJJ     KKKK     LLLL     MMMM     NNNN     OOOO     PPPP     QQQQ     RRRR     SSSS     TTTT     UUUU"
#set list_set     = "AAAAA    BBBBB    CCCCC    DDDDD    EEEEE    FFFFF    GGGGG    HHHHH    IIIII    JJJJJ    KKKKK    LLLLL    MMMMM    NNNNN    OOOOO    PPPPP    QQQQQ    RRRRR    SSSSS    TTTTT    UUUUU"
#set list_set     = "AAAAAA   BBBBBB   CCCCCC   DDDDDD   EEEEEE   FFFFFF   GGGGGG   HHHHHH   IIIIII   JJJJJJ   KKKKKK   LLLLLL   MMMMMM   NNNNNN   OOOOOO   PPPPPP   QQQQQQ   RRRRRR   SSSSSS   TTTTTT   UUUUUU"
#set list_set     = "AAAAAAA  BBBBBBB  CCCCCCC  DDDDDDD  EEEEEEE  FFFFFFF  GGGGGGG  HHHHHHH  IIIIIII  JJJJJJJ  KKKKKKK  LLLLLLL  MMMMMMM  NNNNNNN  OOOOOOO  PPPPPPP  QQQQQQQ  RRRRRRR  SSSSSSS  TTTTTTT  UUUUUUU"
#set list_set     = "AAAAAAAA BBBBBBBB CCCCCCCC DDDDDDDD EEEEEEEE FFFFFFFF GGGGGGGG HHHHHHHH IIIIIIII JJJJJJJJ KKKKKKKK LLLLLLLL MMMMMMMM NNNNNNNN OOOOOOOO PPPPPPPP QQQQQQQQ RRRRRRRR SSSSSSSS TTTTTTTT UUUUUUUU"
#set list_xsid    = "1 2 3 4 5 6 7 8" # for xsspin
set list_xsid    = "1 2 3 4 5 6"
#set list_xsid    = "3 6"

set outdir_merge = `dirname ${3}`/`basename ${3}`_merge/

mkdir -p $3;
mkdir -p ${outdir_merge};
     
foreach set($list_set)
  foreach xsid($list_xsid)
    (./bcs_lr_lep_sig $1 $set $xsid 0 $2 $3 $4     >> log/log_bcs_lr_lep_sig_xsid${xsid}_lep0_${4}.log ) >>& error.log
    (./bcs_lr_lep_sig $1 $set $xsid 1 $2 $3 $4     >> log/log_bcs_lr_lep_sig_xsid${xsid}_lep1_${4}.log ) >>& error.log
  end

  (./bcs_merge  $1 $set $3 ${outdir_merge}     >> log/log_bcs_merge_set${set}.log              ) >>& error.log
end
