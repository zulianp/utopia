#!/bin/bash
set -e
set -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic:.\n The output is logged in make_basic.log.\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2

_basic_build(){

	touch make_basic.log

	cmake $SCRIPT_DIR/../ -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_TRILINOS=OFF -DUTOPIA_ENABLE_PETSC=OFF -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_TESTS=ON -DCMAKE_INSTALL_PREFIX=$prefix | tee make_basic.log

	read -p "Continue ? y/n" -n 1 -r
	echo    # (optional) move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
	    # do dangerous stuff
	    make -j$N_THREADS complete | tee -a make_basic.log
		./utopia_bench | tee -a make_basic.log
		./utopia_test | tee -a make_basic.log
		make install | tee -a make_basic.log
		make -j$N_THREADS test_install | tee -a make_basic.log
	fi
}

if [[ -d build_basic ]]
then
	cd build_basic
	rm -rf *
	_basic_build
else
	if [[ ! -d build_basic ]]
	then
		mkdir build_basic
		cd build_basic
		rm -rf *
		_basic_build
	fi
fi