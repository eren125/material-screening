#!/bin/bash

# Reads all variables defined in set_environment and forwards their definition to $1
THISDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
NAMES=( $(grep -oP "(?<=export )([^=]*)" ${THISDIR}/set_environment) )
VARS="`set -o posix ; set`"
source ${THISDIR}/set_environment
SCRIPT_VARS="`grep -vFe "$VARS" <<<"$(set -o posix ; set)" | grep -v ^VARS=`"
unset VARS
for x in "${NAMES[@]}"
do
    echo `echo $SCRIPT_VARS | grep -oP "\b$x\b=[^ ]*"` >> $1
done
