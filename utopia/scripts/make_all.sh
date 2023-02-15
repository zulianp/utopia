#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Cmake All Script:\n"


_build_all_vc(){
	cmake .. -DUTOPIA_ENABLE_VC
}

# Non VC Environment

_build_all(){
	cmake .. -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_SLEPC=ON -DUTOPIA_ENABLE_METIS=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PARMETIS=OFF -DUTOPIA_ENABLE_TRACE=ON -DUTOPIA_ENABLE_SCRIPTING=OFF -DUTOPIA_PULL_REQUEST_MODE=ON -DUTOPIA_ENABLE_SANITIZER=ON -DTrilinos_DIR=/Users/dylan/Documents/Summer-Internship/Installations/Trilinos/lib/cmake/Trilinos
}

if [[ -d build ]]
then
	cd build
	_build_all
fi

if [[ ! -d build ]]
then
	mkdir build
	cd build
	_build_all
fi