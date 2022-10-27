#!/bin/bash
source ./set_environment
cd "$(dirname "$0")"
ESCAPED_RASPA_DIR=`echo $RASPA_DIR | sed "s=\/=\\\/=g"`
ls $RASPA_DIR/share/raspa/structures/cif/* \
| sed -e '1iStructures' \
| sed -e "s/$ESCAPED_RASPA_DIR\/share\/raspa\/structures\/cif\///g" > structures.csv
