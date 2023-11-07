#!/bin/bash
set -e
set -o pipefail

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Local Install of dependencies petsc and trilinos:\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2

# Basic build, petsc, lapack.
_local_install() {
	touch make_local_install.log

	cmake .. -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TRILINOS=ON -DCMAKE_INSTALL_PREFIX=$prefix | tee make_local_install.log

	read -p "Continue ? y/n" -n 1 -r
	echo "\n"    # (optional) move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		#Make petsc	
		make -j$N_THREADS petsc | tee -a make_local_install.log

	  	#Make Trilinos
		make -j$N_THREADS trilinos | tee -a make_local_install.log

		# Another cmake to find dependencies.
		cmake .. | tee -a make_local_install.log

		make -j$N_THREADS | tee make_local_install.log
		make install | tee make_local_install.log
		make -j$N_THREADS test_install | tee make_local_install.log
	fi
}

if [[ -d build_local_install ]]
then
	cd build_local_install
	rm -rf *
	_local_install
else
	if [[ ! -d build_local_install ]]
	then
		mkdir build_local_install
		cd build_local_install
		rm -rf *
		_local_install
	fi
fi