#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script:\n"


# For now just simple build with BLAS, YAML_CPP, MPI, PETSC, TRILINOS turned on.
# Add more down the line. 
_simple_build(){
	cmake .. -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_BLAS=ON -DTrilinos_DIR=/Users/dylan/Documents/Summer-Internship/Installations/Trilinos/lib/cmake/Trilinos
}

_fluya_mode(){
	cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DTrilinos_DIR=/Users/dylan/Documents/Summer-Internship/Installations/Trilinos/lib/cmake/Trilinos
}

if [[ -d build ]]
then
	cd build
	_simple_build
fi

if [[ ! -d build ]]
then
	mkdir build
	cd build
	_simple_build
fi

