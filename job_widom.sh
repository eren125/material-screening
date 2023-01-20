#!/bin/bash
source ./set_environment
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

$MATSCREEN_PYTHON $MATSCREEN/screen.py -ppn 20 -n 20 -s $MATSCREEN/data/structures_sym.csv -m xenon -g yes -N 100000 -t widom_nogrid

rm -rf Movies/
rm -rf VTK/
rm -rf Restart/
rmdir Scripts
