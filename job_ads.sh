#!/bin/bash
export MATSCREEN=/home/emmanuel/Documents/material-screening
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $MATSCREEN/screen.py -ppn 20 -n 20 -s $MATSCREEN/data/structures_sample.csv -m xenon -N 100 -t ads 

rm -rf Movies/
rm -rf VTK/
rm -rf Restart/
rmdir Scripts
