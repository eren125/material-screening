#!/bin/bash
shopt -s extglob
if [ -f 'DATA' ]; then rm DATA; fi
touch DATA
for line in `ls Output/System_0/output_*`
do
  struc=${line/Output\/System\_0\/output\_}
  struc=${struc/.data}
  grep -m 1 'Host\/Adsorbate energy:' $line | sed -e "s/Host\/Adsorbate energy:/$struc/g" >> DATA
done
