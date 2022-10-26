#!/bin/bash
shopt -s extglob
if [ -f 'DATA' ]; then rm DATA; fi
touch DATA
i=1
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/\_+([0-99]).+([0-99]).+([0-99])\_298.000000_*.data}
  grep 'Average loading absolute \[mol/' $line | sed -e "s/Average loading absolute \[mol\/kg framework\]/$struc/g;s/+\/-//g;s/\[-\]//g" >> DATA 
done
