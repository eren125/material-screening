#!/bin/bash
#TODO integrate it into the python script as an important environment variable
start=$SECONDS
export CURRENTDIR=$PWD
export MATSCREEN=/home/emmanuel/Documents/coudert_lab/material-screening
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $MATSCREEN/screen.py -ppn 1 -n 1 -s $MATSCREEN/data/structures.csv -m xenon -t pore -f UFF -X glost -r 2
echo $(( SECONDS - start ))
