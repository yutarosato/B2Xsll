#!/bin/csh -f

#make;


if( 0 == 1 ) then
  foreach axis (1 3 5 7 9 10) # 1(D0), 3(D*), 5(delta-M), 7(Xd), 9(Xd*), 10(D0-rev)
    bsub -q e ./exe_var.sh  $axis  # [fl_axis]
  end
endif

if( 0 == 1 ) then
  foreach par (0 1) #0(K), 1(pi)
  bsub -q e ./exe_pid.sh 1 $par # [fl_axis(pid)] [fl_par] # mu-id
  bsub -q e ./exe_pid.sh 2 $par # [fl_axis(pid)] [fl_par] #  e-id
  end
endif

if( 0 == 1 ) then
  foreach par (0 1) #0(K), 1(pi)
    bsub -q e ./exe_p_cos.sh 0 $par 0 # [fl_axis] [fl_par] [fl_adef] # mu-id
    bsub -q e ./exe_p_cos.sh 1 $par 0 # [fl_axis] [fl_par] [fl_adef] #  e-id
    bsub -q e ./exe_p_cos.sh 0 $par 1 # [fl_axis] [fl_par] [fl_adef] # mu-id
    bsub -q e ./exe_p_cos.sh 1 $par 1 # [fl_axis] [fl_par] [fl_adef] #  e-id
  end
endif

set adef  = 0 # [01]
if( 1 == 1 ) then
  foreach par (0 1) #0(K), 1(pi)
    foreach mom ( `seq 0 15`) # for adef=0
#    foreach mom ( `seq 0 10`) # for adef=1
      foreach cos ( `seq 0 10`)
#         bsub -q e ./exe_d0_plot.sh 0 $par $mom $cos $adef # [fl_pid] [fl_par] [fl_mom] [fl_cos][fl_adef] // mu-id
#         bsub -q e ./exe_d0_plot.sh 1 $par $mom $cos $adef                                               //  e-id
#         bsub -q e ./exe_d0_plot.sh 2 $par $mom $cos $adef                                               // no-pid
         bsub -q e ./exe_d0_plot_data.sh 0 $par $mom $cos $adef # [fl_pid] [fl_par] [fl_mom] [fl_cos][fl_adef] // mu-id
         bsub -q e ./exe_d0_plot_data.sh 1 $par $mom $cos $adef                                               //  e-id
         bsub -q e ./exe_d0_plot_data.sh 2 $par $mom $cos $adef                                               // no-pid
      end
    end
  end
endif
