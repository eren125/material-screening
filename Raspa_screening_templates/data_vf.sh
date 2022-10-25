#!/bin/bash
shopt -s extglob
if [ -f 'DATA' ]; then rm DATA; fi
touch DATA
i=1
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/\_+([0-9])*.*.*\_298.000000\_*.data}
  grep 'Rosenbluth factor new:' $line | sed -e "s/Rosenbluth factor new:/$struc/g;s/\[-\]//g" >> DATA
done
