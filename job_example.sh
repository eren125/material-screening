#!/bin/bash
#TODO integrate it into the python script as an important environment variable
start=$SECONDS
export CURRENTDIR=$PWD
export MATSCREEN=/home/emmanuel/Documents/material-screening
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $MATSCREEN/screen.py -ppn 20 -n 20 -s $MATSCREEN/data/structures.csv -m xenon -t sample

echo $(( SECONDS - start ))
