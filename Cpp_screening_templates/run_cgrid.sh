#!/bin/bash 
node=$(hostname)
source ${MATSCREEN}/set_environment

cd "$(dirname "$0")" 

STRUC_DIR=$RASPA_DIR/share/raspa/structures/cif

ads_grid.out $STRUC_DIR/$1.cif $RASPA_DIR/share/raspa/forcefield/FORCE_FIELD/force_field_mixing_rules.def $3 CUTOFF ATOMS TIMESTEP 20 0.8 >> PATH/cpp_output_N_cycles.csv
