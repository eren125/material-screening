#!/bin/bash
cd "$(dirname "$0")"
ls ~/Simulations/share/raspa/structures/cif/* \
| sed -e '1iStructures' \
| sed -e "s/\/home\/emmanuel\/Simulations\/share\/raspa\/structures\/cif\///g" > structures.csv
