#!/usr/bin/env bash

# WARNING! If Trilinos complains about getopt you need to make sure that 
# the update version is in the path!

set -e

export PATH=$PWD:$PATH

if [[ -z "$UTOPIA_FE_EXEC" ]]
then
	echo "Error! Please define UTOPIA_FE_EXEC=<path_to_utopia_exectuable>"
	exit -1
fi

if [[ -z "$TRILINOS_DIR" ]]
then
	echo "Error! Please define TRILINOS_DIR=<path_to_trilinos_installation>"
	exit -1
fi

if (($# != 3))
then
	printf "usage: $0 <fluid.e> <solid.e> <n paritions>\n" 1>&2
	exit -1
fi

FLUID_MESH=fluid.e
SOLID_MESH=solid.e

# Copy the mesh to local folder
cp $1 $FLUID_MESH
cp $2 $SOLID_MESH

N_PARTS=$3
RESCALE_IMBALANCE=1

rm -f resample.yaml temp_.yaml
( echo "cat <<EOF >resample.yaml";
  cat static_load_balancing_tpl.yaml;
  echo "EOF";
) >temp_.yaml
. temp_.yaml

rm temp_.yaml
cat resample.yaml

# Parallel node_to_element_matrix does not work yet
# $TRILINOS_DIR/bin/decomp -p $N_PARTS $FLUID_MESH
# $TRILINOS_DIR/bin/decomp -p $N_PARTS $SOLID_MESH

# mpiexec -np $N_PARTS \
$UTOPIA_FE_EXEC @file resample.yaml

# costs=(`ls cost.*.*`)
# $TRILINOS_DIR/bin/epu -auto ${costs[0]}

rm -f input-ldbl temp_.txt
( echo "cat <<EOF >input-ldbl";
  cat input-ldbl-tpl.txt;
  echo "EOF";
) >temp_.txt
. temp_.txt

rm temp_.txt

$TRILINOS_DIR/bin/decomp -p $N_PARTS -V -i input-ldbl $FLUID_MESH


# clean-up
rm resample.yaml
rm input-ldbl

