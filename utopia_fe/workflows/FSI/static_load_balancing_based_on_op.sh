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

if [[ -z "$SFEM_DIR" ]]
then
	echo "Error! Please define SFEM_DIR=<path_to_trilinos_installation>"
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

EXPORT_EXAMPLE_COUPLED_SYSTEM=true

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

raw_system=raw_system
raw_decomp=decomp.raw
raw_mesh=raw_mesh
text_decomp=decomp.txt

set -x
$SFEM_DIR/python/mesh/exodusII_to_raw.py $FLUID_MESH $raw_mesh
$SFEM_DIR/workflows/decomposition/partition_mesh_based_on_operator.sh $raw_mesh $raw_system $N_PARTS $raw_decomp
$SFEM_DIR/python/raw2text.py $raw_decomp int32 $text_decomp

$TRILINOS_DIR/bin/decomp -p $N_PARTS -l infile=$text_decomp $FLUID_MESH


$SFEM_DIR/python/fp_convert.py  $raw_decomp proc-float32.raw int32 float32
$SFEM_DIR/python/mesh/raw_to_db.py $raw_mesh decomp-vis.vtk --cell_data="proc-float32.raw" --cell_data_type="float32"

# clean-up
rm resample.yaml


