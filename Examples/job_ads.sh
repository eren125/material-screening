#!/bin/bash
source ${MATSCREEN}/set_environment
#export NODES=$(srun hostname | sort | uniq)
#echo "Running on nodes:"
#echo "${NODES}"

$MATSCREEN_PYTHON $MATSCREEN/screen.py -ppn 20 -n 20 -s $MATSCREEN/data/structures_sample.csv -th 20 -m xenon -N 100 -t ads

rm -rf Movies/
rm -rf VTK/
rm -rf Restart/
rmdir Scripts
