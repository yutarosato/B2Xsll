#!/bin/csh -f

(./draws_2d_likelihood      ${1} ${2} ${3} ${4} ${5} ${6} ${7} 0 > log/log_2d_${7}_${6}_s0${1}_set${2}.log      ) >>& error.log
#(./draws_2d_likelihood_proj ${1} ${2} ${3} ${4} ${5} ${6} ${7} 0 > log/log_2d_${7}_${6}_proj_s0${1}_set${2}.log ) >>& error.log
