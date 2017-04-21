#!/bin/csh -f

#set DIR = "20121117_BGM/fitter_test/toy/nsig/"
set DIR = "log/"


echo "[CHECK FITTING STATUS]"
grep STATUS     ${DIR}/log_Mbc_bin2_toy_* | grep -v INITIATE | grep -v OK | grep -v CONVERGED
grep -i limit   ${DIR}/log_Mbc_bin2_toy_*
#grep -i warning ${DIR}/log_Mbc_bin2_toy_* | grep -v -i Plotting;
#grep -i error   ${DIR}/log_Mbc_bin2_toy_* | grep -v  "ERROR MATRIX" | grep -v SIZE | grep -v "pdf_peak_swap" | grep -v "pdf_peak_cc";

echo "[MAKE LINEARITY-FILES]"
set kakunin = $<

echo "[tmp0_nq2.dat ]"; cat ${DIR}/log_Mbc_bin2_toy_* | grep hooge0  | tee tmp0_nq2.dat;  echo "*" >> tmp0_nq2.dat;
echo "[tmp1_nq2.dat ]"; cat ${DIR}/log_Mbc_bin2_toy_* | grep hooge1  | tee tmp1_nq2.dat;  echo "*" >> tmp1_nq2.dat;
echo "[tmp2_afb.dat ]"; cat ${DIR}/log_Mbc_bin2_toy_* | grep hoooge2 | tee tmp2_afb.dat;  echo "*" >> tmp2_afb.dat;
awk '{if(NR%4==0){tot+=$3; totE+=$4*$4; print $1" "tot" "sqrt(totE);tot=0;totE=0;}else{ tot+=$3; totE+=$4*$4;}}' tmp0_nq2.dat | tee tmp0_ntot.dat; echo "*" >> tmp0_ntot.dat;
awk '{if(NR%4==0){tot+=$3; totE+=$4*$4; print $1" "tot" "sqrt(totE);tot=0;totE=0;}else{ tot+=$3; totE+=$4*$4;}}' tmp1_nq2.dat | tee tmp1_ntot.dat; echo "*" >> tmp1_ntot.dat;

echo "[SUBMIT]"
set kakunin = $<

bsub -q e ./exe_linearity_check_bin2_toy.sh 0;
bsub -q e ./exe_linearity_check_bin2_toy.sh 1;
