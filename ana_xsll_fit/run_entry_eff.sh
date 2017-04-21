#!/bin/csh -f

make draws_1d_entry_eff;
cat log_gen/log_entry_gen_*        | grep HOGE > tmp_gen.dat; echo "*" >> tmp_gen.dat;
#cat log_rec_mbc527/log_entry_rec_* | grep HOGE > tmp_rec.dat; echo "*" >> tmp_gen.dat;
cat log_rec_mbc522/log_entry_rec_* | grep HOGE > tmp_rec.dat; echo "*" >> tmp_gen.dat;

./draws_1d_entry_eff;

