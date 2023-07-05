#!/usr/bin/env bash

# WARNING! If Trilinos complains about getopt you need to make sure that 
# the update version is in the path!

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

export PATH=$SCRIPTPATH:$PATH

if [[ -z "$TRILINOS_DIR" ]]
then
	echo "Error! Please define TRILINOS_DIR=<path_to_trilinos_installation>"
	exit -1
fi


if (($# != 3))
then
	printf "usage: $0 <mesh.e> <n paritions> <decomp.txt>\n" 1>&2
	exit -1
fi

workspace=workspace
mkdir -p $workspace

MESH=$workspace/"`basename $1`"

set -x

# Copy the mesh to local folder
cp $1 $MESH
N_PARTS=$2
DECOMP_FILE=$3

$TRILINOS_DIR/bin/decomp -p $N_PARTS -l infile=${DECOMP_FILE} $MESH
