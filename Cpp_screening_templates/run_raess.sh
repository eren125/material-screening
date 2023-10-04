#!/bin/bash
source ${MATSCREEN}/set_environment

cd "$(dirname "$0")"

STRUC_DIR=$RASPA_DIR/share/raspa/structures/cif

raess $STRUC_DIR/$1.cif $RASPA_DIR/share/raspa/forcefield/FORCE_FIELD/force_field_mixing_rules.def $3 CUTOFF N_cycles ATOMS REJECT TIMESTEP >> PATH/cpp_output_N_cycles_REJECT_TIMESTEP.csv
