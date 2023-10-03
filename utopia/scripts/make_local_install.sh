#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Local Install of dependencies petsc and trilinos:\n"


# Basic build, petsc, lapack.
_local_install(){
	touch make_local_install.log

	cmake .. -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TRILINOS=ON -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia_local_install | tee make_local_install.log


	#Make petsc	
	make -j$1 petsc | tee -a make_local_install.log

	# Define ENV variables for Petsc. Should be defined in .bashrc.
	export PETSC_DIR="$PWD/../../external/petsc"
  	export PETSC_ARCH=arch-darwin-c-debug 
  	export SLEPC_DIR="$PWD/../../external/petsc"

  	#Make Trilinos
	make -j$1 trilinos | tee -a make_local_install.log

	# Another cmake to find dependencies.
	cmake .. | tee -a make_local_install.log

	make -j$1 | tee make_local_install.log
	make install | tee make_local_install.log
	make -j$1 test_install | tee make_local_install.log

}

jobs=$@

if [[ -d build_local_install ]]
then
	cd build_local_install
	rm -rf *
	_local_install $jobs
fi

if [[ ! -d build_local_install ]]
then
	mkdir build_local_install
	cd build_local_install
	rm -rf *
	_local_install $jobs
fi