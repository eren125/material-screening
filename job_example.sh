#!/bin/bash
export RASPASCREEN=/home/emmanuel/Documents/These_1/08_Raspa_screening 
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $RASPASCREEN/raspascreen.py -p 20 -n 20 -s $RASPASCREEN/data/structures_done.csv -m xenon -ps 0.1 1 10 1e2 1e3 1e4 1e5 -N 100000 -t ads -o .
