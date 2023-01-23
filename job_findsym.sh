#!/bin/bash
start=$SECONDS
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

python3 $MATSCREEN/screen.py -ppn 1 -n 1 -s $MATSCREEN/data/structures.csv -t findsym -X glost
echo $(( SECONDS - start ))
