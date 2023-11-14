#!/bin/bash
set -e
set -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script build all:\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2

# Non VC Environment

_build_all(){
	touch make_all.log

	cmake $SCRIPT_DIR/../ -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_BENCHMARK=ON -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_SLEPC=ON -DUTOPIA_ENABLE_METIS=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PARMETIS=OFF -DUTOPIA_ENABLE_TRACE_EXPR=OFF -DUTOPIA_ENABLE_SCRIPTING=OFF -DUTOPIA_PULL_REQUEST_MODE=ON -DCMAKE_INSTALL_PREFIX=$prefix | tee make_all.log
	
	read -p "Continue ? y/n" -n 1 -r
	echo    # (optional) move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		make -j$N_THREADS complete | tee make_all.log
		make install | tee make_all.log
		./utopia_bench | tee make_all.log
		./utopia_test | tee make_all.log
		make -j$N_THREADS test_install | tee make_all.log
	fi
}

if [[ -d build_all ]]
	then
		cd build_all
		_build_all
	else
	if [[ ! -d build_all ]]
	then
		mkdir build_all
		cd build_all
		_build_all
	fi
fi
