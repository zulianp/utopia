#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic:\n"


# Basic build, petsc, lapack.
_basic_build(){
	rm CMakeCache.txt
	cmake .. -DUTOPIA_ENABLE_TRILINOS=OFF -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_YAML_CPP=OFF -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PETSC=OFF -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_BENCHMARK=ON
	make -j 
	./utopia_bench
	./utopia_test
}

if [[ -d build_basic ]]
then
	cd build_basic
	_basic_build
fi

if [[ ! -d build_basic ]]
then
	mkdir build_basic
	cd build_basic
	_basic_build
fi