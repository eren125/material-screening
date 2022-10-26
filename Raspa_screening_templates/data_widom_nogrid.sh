#!/bin/bash
shopt -s extglob
if [ -f 'Henry.txt' ]; then rm Henry.txt; fi
if [ -f 'Energy.txt' ]; then rm Energy.txt; fi
touch Henry.txt
touch Energy.txt
i=1
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/\_+([0-99]).+([0-99]).+([0-99])\_298.000000_0.data}
  grep "\[xenon\] Average Henry coefficient:" $line | sed -e "s/\[xenon\] Average Henry coefficient:/$struc/g" >> Henry.txt
  grep "\[xenon\] Average  <U_gh>_1-<U_h>_0:" $line | sed -e "s/\[xenon\] Average  <U_gh>_1-<U_h>_0:/$struc/g" >> Energy.txt
  grep "total time:" $line | sed -e "s/total time:/$struc/g" >> time.txt
done

