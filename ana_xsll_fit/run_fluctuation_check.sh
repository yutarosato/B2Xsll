#!/bin/csh -f

make fluctuation_check;

#set DIR = "log/"
#set DIR = "20130101/results_syst_sigshape/"
#set DIR = "20130101/results_syst_cf_width/"
set DIR = "20130101/results_syst_scf/"

echo "[MAKE FLUCTUATION-FILES]"
#set kakunin = $<
echo "[tmp0.dat]"; cat ${DIR}/log_Mbc_*_0q2* | grep hoooge2 | tee tmp0.dat
echo "[tmp1.dat]"; cat ${DIR}/log_Mbc_*_1q2* | grep hoooge2 | tee tmp1.dat
echo "[tmp2.dat]"; cat ${DIR}/log_Mbc_*_2q2* | grep hoooge2 | tee tmp2.dat
echo "[tmp3.dat]"; cat ${DIR}/log_Mbc_*_3q2* | grep hoooge2 | tee tmp3.dat

echo "[SUBMIT]"
#set kakunin = $<

./fluctuation_check | tee log_fluctuation.log
echo $DIR >> log_fluctuation.log
