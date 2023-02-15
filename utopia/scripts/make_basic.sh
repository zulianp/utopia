#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic:\n"


# Basic build, petsc, lapack.
_basic_build(){
	rm CMakeCache.txt
	cmake .. -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON
}

if [[ -d build ]]
then
	cd build
	_basic_build
fi

if [[ ! -d build ]]
then
	mkdir build
	cd build
	_basic_build
fi