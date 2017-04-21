#!/bin/csh -f

foreach exp(`cat ~/script/nBB2.txt`)
foreach mc(0 1 2 3 4)
#foreach mc(0)
   bsub -q ex ./exe_entry_rec.sh $mc $exp 0 0
end
end
