#!/bin/bash
shopt -s extglob
if [ -f 'Loading.txt' ]; then rm Loading.txt; fi
if [ -f 'Enthalpy_temp.txt' ]; then rm Enthalpy_temp.txt; fi
touch Loading.txt
touch Enthalpy_temp.txt
i=1
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/\_+([0-9])*.+([0-9])*.+([0-9])*\_298.000000\_*.data}
  grep 'Average loading absolute \[mol/' $line | sed -e "s/Average loading absolute \[mol\/kg framework\]/$struc/g;s/+\/-//g;s/\[-\]//g" >> Loading.txt 
  echo $struc >> Enthalpy_temp.txt
  grep -A 64 'Enthalpy' $line >> Enthalpy_temp.txt 
done
sed -n '1~66p;14~66p;27~66p;40~66p' Enthalpy_temp.txt > Enthalpy.txt
rm Enthalpy_temp.txt
