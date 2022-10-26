#!/bin/bash
#TODO integrate it into the python script as an important environment variable
start=$SECONDS
MATSCREEN=/home/emmanuel/Documents/material-screening # replace with your path to material-screening
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

$MATSCREEN/screen.py -ppn 20 -n 20 -s $MATSCREEN/data/structures.csv -m xenon -t sample

echo $(( SECONDS - start ))
