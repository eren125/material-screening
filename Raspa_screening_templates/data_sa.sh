#!/bin/bash
shopt -s extglob
if [ -f 'sa.txt' ]; then rm sa.txt; fi
touch sa.txt
i=1
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/\_+([0-99]).+([0-99]).+([0-99])\_298.000000_0.data}
  grep -A 2 "Average surface area:" $line | sed -e "s/Average surface area:/$struc/g" | tr '\n' ' ' >> sa.txt
  echo -e "" >> sa.txt
done
