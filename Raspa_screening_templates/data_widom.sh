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
  struc=${struc/\_+([0-9])*.+([0-9])*.+([0-9])*\_298.000000*.data}
  grep "\[MOLECULE\] Average Henry coefficient:" $line | sed -e "s/\[MOLECULE\] Average Henry coefficient:/$struc/g" >> Henry.txt
  grep "\[MOLECULE\] Average  <U_gh>_1-<U_h>_0:" $line | sed -e "s/\[MOLECULE\] Average  <U_gh>_1-<U_h>_0:/$struc/g" >> Energy.txt
done
