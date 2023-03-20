#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic:\n"


# Basic build, petsc, lapack.
_basic_build(){
	cmake .. -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_BENCHMARK=ON -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_SCRIPTING=ON -DUTOPIA_ENABLE_SLEPC=ON -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia

	make -j complete
	./utopia_bench
	./utopia_test
	make install
}

if [[ -d build_basic ]]
then
	cd build_basic
	rm -rf *
	_basic_build
fi

if [[ ! -d build_basic ]]
then
	mkdir build_basic
	cd build_basic
	rm -rf *
	_basic_build
fi