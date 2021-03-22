#!/bin/bash
export RASPASCREEN=/home/emmanuel/Documents/raspa-screener
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $RASPASCREEN/raspascreen.py -ppn 20 -n 20 -s $RASPASCREEN/data/structures_sample.csv -m xenon krypton -N 100000 -t grid 
