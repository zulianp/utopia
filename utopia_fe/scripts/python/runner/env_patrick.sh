#!/bin/bash

export CODE_DIR=/Users/patrickzulian/Desktop/code
export UTOPIA_FE_SRC_DIR=$CODE_DIR/utopia/utopia_fe/
export UTOPIA_FE_EXEC_PYTHON=$UTOPIA_FE_SRC_DIR/scripts/python/runner
export UTOPIA_FE_EXEC=$UTOPIA_FE_SRC_DIR/build_petsc/utopia_fe_exec

export TRILINOS_DECOMP=$CODE_DIR/installations/Trilinos/bin/decomp
export TRILINOS_EPU=$CODE_DIR/installations/Trilinos/bin/epu
export UTOPIA_LAUNCH_CMD=mpiexec