#!/bin/bash
#TODO integrate it into the python script as an important environment variable
start=$SECONDS
export CURRENTDIR=$PWD
export MATSCREEN=/home/emmanuel/Documents/coudert_lab/material-screening
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $MATSCREEN/screen.py -ppn 1 -n 1 -s $MATSCREEN/data/structures_sym.csv -m xenon -t raess -f UFF -g yes -N 2000 -rj 0.85 -r 1.6 -o .
bash $MATSCREEN/copy_env.sh set_environment
echo $(( SECONDS - start ))
