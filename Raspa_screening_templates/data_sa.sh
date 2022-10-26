#!/bin/bash
shopt -s extglob
if [ -f 'sa.txt' ]; then rm sa.txt; fi
touch sa.txt
i=1
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/\_+([0-99]).+([0-99]).+([0-99])\_298.000000_0.data}
  grep "Surface area:" $line | sed -n 3p | sed -e "s/Surface area:/$struc/g" >> sa.txt
done
