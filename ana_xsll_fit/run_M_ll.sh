#!/bin/csh -f

foreach cut(`seq 0 30`)
  bsub -q ex ./exe_M_ll.sh 0 $cut
  bsub -q ex ./exe_M_ll.sh 1 $cut
end
