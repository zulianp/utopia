#!/usr/bin/env bash

# WARNING! If Trilinos complains about getopt you need to make sure that 
# the update version is in the path!

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

export PATH=$SCRIPTPATH:$PATH

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

if [[ -z "$LAUNCH" ]]
then
	LAUNCH=mpiexec
fi

if (($# != 3))
then
	printf "usage: $0 <fluid.e> <solid.e> <n paritions>\n" 1>&2
	exit -1
fi

workspace=workspace
mkdir -p $workspace

FLUID_MESH=$workspace/`basename $1`
SOLID_MESH=$workspace/`basename $2`

# Copy the mesh to local folder
cp $1 $FLUID_MESH
cp $2 $SOLID_MESH

N_PARTS=$3
# Control scale of the weights (effective range?)

if [[ -z "$RESCALE_IMBALANCE" ]]
then
	RESCALE_IMBALANCE=2
fi

if [[ -z "$DISPLACE_SOLID" ]]
then
	DISPLACE_SOLID=true
fi

rm -f resample.yaml temp_.yaml
( echo "cat <<EOF >resample.yaml";
  cat $SCRIPTPATH/static_load_balancing_tpl.yaml;
  echo "EOF";
) >temp_.yaml
. temp_.yaml

rm temp_.yaml
cat resample.yaml

# Parallel cost estimation or not
if [[ -z "$COST_ESTIMATION_N_PROCS" ]]
then
	$UTOPIA_FE_EXEC @file resample.yaml
else
	# Parallel node_to_element_matrix does not work yet
	cd $workspace
	$TRILINOS_DIR/bin/decomp -p $COST_ESTIMATION_N_PROCS --spectral -V `basename $FLUID_MESH`
	$TRILINOS_DIR/bin/decomp -p $COST_ESTIMATION_N_PROCS --spectral -V `basename $SOLID_MESH`
	cd -

	$LAUNCH -np $COST_ESTIMATION_N_PROCS \
	$UTOPIA_FE_EXEC @file resample.yaml
	costs=(`ls cost.*.*`)
	$TRILINOS_DIR/bin/epu -auto ${costs[0]}
fi

rm -f input-ldbl temp_.txt
( echo "cat <<EOF >input-ldbl";
  cat $SCRIPTPATH/input-ldbl-tpl.txt;
  echo "EOF";
) >temp_.txt
. temp_.txt

rm temp_.txt

$TRILINOS_DIR/bin/decomp -p $N_PARTS --spectral -V -i input-ldbl $FLUID_MESH


# clean-up
rm resample.yaml
rm input-ldbl

