#!/bin/csh -f

#make;

# if chagne it, should change the infile[012] in draws_1d_d0_fit.cpp
#set DIR = 'pic_d0_plot_xdst050'
#set DIR = 'pic_d0_plot_xdst035'
#set DIR = 'pic_d0_plot_data_xdst050'
set DIR = 'pic_d0_plot_data_xdst035'

set CNT=0

#foreach f($DIR"/d0_plot_data_pid2_par1_mom1_cos1_"*".root") # for debug
foreach f($DIR"/d0_plot_data_pid2_"*".root")
   echo -n $CNT ' '$f
   set file0 = `echo $f | sed -e "s/_pid2_/_pid0_/g"`
   set file1 = `echo $f | sed -e "s/_pid2_/_pid1_/g"`

     if !(-e $file0) then
     echo " : NOT FOUND pid0-file -> Skip"
     continue
     endif
     if !(-e $file1) then
     echo " : NOT FOUND pid1-file -> Skip"
     continue
     endif
     
   @ CNT += 1

   set BASE = `basename $f .root`
####################### MC #######################
#   set PAR  = `echo $BASE | cut -d "_" -f 4 | sed -e "s/par//g"`
#   set MOM  = `echo $BASE | cut -d "_" -f 5 | sed -e "s/mom//g"`
#   set COS  = `echo $BASE | cut -d "_" -f 6 | sed -e "s/cos//g"`
#   set ADEF = `echo $BASE | cut -d "_" -f 7 | sed -e "s/adef//g"`
#   echo -n " : $PAR $MOM $COS $ADEF     -> "
#   bsub -q e ./exe_d0_fit.sh $PAR $MOM $COS $ADEF

###################### DATA ######################
   set PAR  = `echo $BASE | cut -d "_" -f 5 | sed -e "s/par//g"`
   set MOM  = `echo $BASE | cut -d "_" -f 6 | sed -e "s/mom//g"`
   set COS  = `echo $BASE | cut -d "_" -f 7 | sed -e "s/cos//g"`
   set ADEF = `echo $BASE | cut -d "_" -f 8 | sed -e "s/adef//g"`
   echo -n " : $PAR $MOM $COS $ADEF     -> "
   bsub -q e ./exe_d0_fit_data.sh $PAR $MOM $COS $ADEF
end
