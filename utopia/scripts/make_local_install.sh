#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Local Install of dependencies petsc and trilinos:\n"


# Basic build, petsc, lapack.
_basic_build(){
	cmake .. -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=ON -DUTOPIA_ENABLE_PETSC=OFF -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia

}

if [[ -d build_local_install ]]
then
	cd build_local_install
	rm -rf *
	_basic_build
fi

if [[ ! -d build_local_install ]]
then
	mkdir build_local_install
	cd build_local_install
	rm -rf *
	_basic_build
fi