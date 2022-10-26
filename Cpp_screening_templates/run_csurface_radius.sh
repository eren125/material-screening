#!/bin/bash 
node=$(hostname)
export MATSCREEN=$CCCSCRATCHDIR/material-screening
source ${MATSCREEN}/set_environment

cd "$(dirname "$0")" 

STRUC_DIR=$RASPA_DIR/share/raspa/structures/cif

${CPP_DIR}/bin/surface_radius.out $STRUC_DIR/$1.cif $RASPA_DIR/share/raspa/forcefield/FORCE_FIELD/force_field_mixing_rules.def TEMPERATURE CUTOFF N_cycles ATOMS TIMESTEP >> PATH/cpp_output_N_cycles.csv
