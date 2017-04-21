#!/bin/csh -f

(./draws_1d_likelihood_lep5 ${1} ${2} ${3} ${4} ${5} ${6} 0 > log/log_1d_lep5_${6}_s0${1}_set${2}.log ) >>& error.log
