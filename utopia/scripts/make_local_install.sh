#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Local Install of dependencies petsc and trilinos:\n"


# Basic build, petsc, lapack.
_local_install(){
	touch make_local_install.log

	cmake .. -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TRILINOS=ON -DCMAKE_INSTALL_PREFIX=$2 | tee make_local_install.log


	#Make petsc	
	make -j$1 petsc | tee -a make_local_install.log

  	#Make Trilinos
	make -j$1 trilinos | tee -a make_local_install.log

	# Another cmake to find dependencies.
	cmake .. | tee -a make_local_install.log

	make -j$1 | tee make_local_install.log
	make install | tee make_local_install.log
	make -j$1 test_install | tee make_local_install.log

}

jobs=$1
prefix=$2

if [[ -d build_local_install ]]
then
	cd build_local_install
	rm -rf *
	_local_install $jobs $prefix
fi

if [[ ! -d build_local_install ]]
then
	mkdir build_local_install
	cd build_local_install
	rm -rf *
	_local_install $jobs $prefix
fi