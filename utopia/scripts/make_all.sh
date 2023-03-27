#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Cmake All Script:\n"


_build_all_vc(){
	cmake .. -DUTOPIA_ENABLE_VC
}

# Non VC Environment

_build_all(){
	cmake .. -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_BENCHMARK=ON -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_MPI=ON -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_SLEPC=ON -DUTOPIA_ENABLE_METIS=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PARMETIS=OFF -DUTOPIA_ENABLE_TRACE=OFF -DUTOPIA_ENABLE_SCRIPTING=OFF -DUTOPIA_PULL_REQUEST_MODE=ON -DUTOPIA_ENABLE_SANITIZER=ON -DTrilinos_DIR=/Users/dylan/Documents/Summer-Internship/Installations/Trilinos/lib/cmake/Trilinos

	make -j4 complete
	# ./utopia_bench
	# ./utopia_test
}

if [[ -d build_all ]]
then
	cd build_all
	_build_all
fi

if [[ ! -d build_all ]]
then
	mkdir build_all
	cd build_all
	_build_all
fi