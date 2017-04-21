#!/bin/csh -f

make draws_fake_rate

# axis definition "adef" should be changed in draws_fake_rate.cpp
#cat log_d0_fit_xdst035/*_adef0.dat > test.dat ;
#cat log_d0_fit_xdst035/*_adef1.dat > test.dat ;
#cat log_d0_fit_xdst050/*_adef0.dat > test.dat ;
#cat log_d0_fit_xdst050/*_adef1.dat > test.dat ;

#cat log_d0_fit_data_xdst035/*_adef0.dat > test.dat ;
#cat log_d0_fit_data_xdst035/*_adef1.dat > test.dat ;
#cat log_d0_fit_data_xdst050/*_adef0.dat > test.dat ;
#cat log_d0_fit_data_xdst050/*_adef1.dat > test.dat ;


./draws_fake_rate
