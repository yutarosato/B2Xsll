#!/bin/csh -f

#make;


#bsub -q s ./exe_Mbc_cc_fit.sh 0-5


if( 0 == 1 ) then
  foreach func (50 51 52 53 15 151 152 153) # 15 151 50 51, 52 53 152 153
      bsub -q e ./exe_Mbc_nopeak_fit.sh   0-5 6 $func 
  end
endif

if( 0 == 1 ) then
  foreach func (15 151 50 51)
    bsub -q s ./exe_Mbc_emu_fit.sh  0-5 $func 
  end
endif

  
if( 1 == 1 ) then
  foreach func (50) # 15 151 50 51
      bsub -q s ./exe_Mbc_shape.sh   0-5 6 $func 
      #./exe_Mbc_shape.sh   0-5 6 $func 
  end
endif
