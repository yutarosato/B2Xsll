#!/bin/csh -f

foreach exp(`cat ~/script/nBB2.txt`)
foreach mc(0 3 4)
   bsub -q ex ./exe_entry_gen.sh $mc $exp 0 0
end
end

set list = "~/script/caseB-mixedjpsi.list"
#echo $list
foreach exp(`cat ~/script/nBB2.txt`)
foreach mc(1) # mixedjpsi
   set exp2 = `echo $exp | sed 's/^0//'`
   set LINEi = `cat $list | grep -n ".e$exp2" | head -1 | awk -F ":" '{print $1}'`
   set LINEf = `cat $list | grep -n ".e$exp2" | tail -1 | awk -F ":" '{print $1}'`
   set NLINE = `expr $LINEf - $LINEi + 1`
   #echo $exp $exp2 $LINEi $LINEf $NLINE
   bsub -q ex ./exe_entry_gen.sh $mc $exp $NLINE $LINEi
end
end

set list = "~/script/caseB-chargedjpsi.list"
#echo $list
foreach exp(`cat ~/script/nBB2.txt`)
foreach mc(2) # chargedjpsi
   set exp2 = `echo $exp | sed 's/^0//'`
   set LINEi = `cat $list | grep -n ".e$exp2" | head -1 | awk -F ":" '{print $1}'`
   set LINEf = `cat $list | grep -n ".e$exp2" | tail -1 | awk -F ":" '{print $1}'`
   set NLINE = `expr $LINEf - $LINEi + 1`
   #echo $exp $exp2 $LINEi $LINEf $NLINE
   bsub -q ex ./exe_entry_gen.sh $mc $exp $NLINE $LINEi
end
end
