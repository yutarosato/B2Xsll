#!/bin/csh -f

make linearity_check;

#set DIR = "trial_8b/"
set DIR = "log/"


echo "[CHECK FITTING STATUS]"
grep STATUS ${DIR}/log_Mbc_func* | grep -v INITIATE | grep -v OK | grep -v CONVERGED
grep -i limit   ${DIR}/log_Mbc_func*
grep -i warning ${DIR}/log_Mbc_func* | grep -v -i Plotting;
grep -i error   ${DIR}/log_Mbc_func* | grep -v  "ERROR MATRIX" | grep -v SIZE;

echo "[MAKE LINEARITY-FILES]"
set kakunin = $<
echo "[tmp0.dat]"; cat ${DIR}/log_Mbc_func* | grep HOGE0 | tee tmp0.dat
echo "[tmp1.dat]"; cat ${DIR}/log_Mbc_func* | grep HOGE1 | tee tmp1.dat
echo "[tmp2.dat]"; cat ${DIR}/log_Mbc_func* | grep HOGE2 | tee tmp2.dat

echo "[SUBMIT]"
set kakunin = $<

bsub -q s ./exe_linearity_check.sh 0;
bsub -q s ./exe_linearity_check.sh 1;
bsub -q s ./exe_linearity_check.sh 2;
