#!/bin/csh -f

make linearity_check;

#set DIR = "test/"
set DIR = "log/"


echo "[CHECK FITTING STATUS]"
grep STATUS ${DIR}/log_Mbc_bin_func* | grep -v INITIATE | grep -v OK | grep -v CONVERGED
grep -i limit   ${DIR}/log_Mbc_bin_func*
grep -i warning ${DIR}/log_Mbc_bin_func* | grep -v -i Plotting;
grep -i error   ${DIR}/log_Mbc_bin_func* | grep -v  "ERROR MATRIX" | grep -v SIZE;

echo "[MAKE LINEARITY-FILES]"
set kakunin = $<
echo "[tmp0_ntot.dat]"; cat ${DIR}/log_Mbc_bin_func* | grep HOGE0   | tee tmp0_ntot.dat; echo "*" >> tmp0_ntot.dat;
echo "[tmp1_ntot.dat]"; cat ${DIR}/log_Mbc_bin_func* | grep HOGE1   | tee tmp1_ntot.dat; echo "*" >> tmp1_ntot.dat;
echo "[tmp2_ntot.dat]"; cat ${DIR}/log_Mbc_bin_func* | grep HOGE2   | tee tmp2_ntot.dat; echo "*" >> tmp2_ntot.dat;
echo "[tmp0_nq2.dat ]"; cat ${DIR}/log_Mbc_bin_func* | grep HOOGE0  | tee tmp0_nq2.dat;  echo "*" >> tmp0_nq2.dat;
echo "[tmp1_nq2.dat ]"; cat ${DIR}/log_Mbc_bin_func* | grep HOOGE1  | tee tmp1_nq2.dat;  echo "*" >> tmp1_nq2.dat;
echo "[tmp2_nq2.dat ]"; cat ${DIR}/log_Mbc_bin_func* | grep HOOGE2  | tee tmp2_nq2.dat;  echo "*" >> tmp2_nq2.dat;
echo "[tmp0_afb.dat ]"; cat ${DIR}/log_Mbc_bin_func* | grep HOOOGE0 | tee tmp0_afb.dat;  echo "*" >> tmp0_afb.dat;
echo "[tmp1_afb.dat ]"; cat ${DIR}/log_Mbc_bin_func* | grep HOOOGE1 | tee tmp1_afb.dat;  echo "*" >> tmp1_afb.dat;
echo "[tmp2_afb.dat ]"; cat ${DIR}/log_Mbc_bin_func* | grep HOOOGE2 | tee tmp2_afb.dat;  echo "*" >> tmp2_afb.dat;

echo "[SUBMIT]"
set kakunin = $<

bsub -q e ./exe_linearity_check_bin.sh 0;
bsub -q e ./exe_linearity_check_bin.sh 1;
bsub -q e ./exe_linearity_check_bin.sh 2;
