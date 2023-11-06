#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script build all:\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

# Non VC Environment
# NEED TO TEST PARMETIS STILL.
# KSP_CONVERGED_CG_NEG_CURVE' is deprecated.
# Sanitizer gives linker error.

_build_all(){

	echo $1

	touch make_all.log
	cmake .. -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_BENCHMARK=ON -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_SLEPC=ON -DUTOPIA_ENABLE_METIS=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PARMETIS=OFF -DUTOPIA_ENABLE_TRACE_EXPR=OFF -DUTOPIA_ENABLE_SCRIPTING=OFF -DUTOPIA_PULL_REQUEST_MODE=ON -DCMAKE_INSTALL_PREFIX=$2 | tee make_all.log

	make -j$1 complete | tee make_all.log
	make install | tee make_all.log
	./utopia_bench | tee make_all.log
	./utopia_test | tee make_all.log
	make -j$1 test_install | tee make_all.log
}

jobs=$N_THREADS
prefix=$2

if [[ -d build_all ]]
then
	cd build_all
	rm -rf *
	_build_all $jobs $prefix
fi

if [[ ! -d build_all ]]
then
	mkdir build_all
	cd build_all
	rm -rf *
	_build_all $jobs $prefix
fi