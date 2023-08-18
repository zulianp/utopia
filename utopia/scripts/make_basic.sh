#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic: Download and install Petsc, Download and instal Trilinos, make complete, benchmark, test, install and finally test_install.\n The output is logged in make_basic.log.\n"

_basic_build(){

	touch make_basic.log

	cmake .. -DUTOPIA_ENABLE_TRACE_EXPR=ON -DUTOPIA_ENABLE_TRACE=ON -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_YAML_CPP=ON -DUTOPIA_ENABLE_POLYMORPHIC=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_BENCHMARK=ON -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_SCRIPTING=OFF -DUTOPIA_ENABLE_SLEPC=ON -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia -DUTOPIA_ENABLE_ISOLVER=ON -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=ON | tee make_basic.log

	
	#Make petsc	
	make -j4 petsc | tee -a make_basic.log

	# Define ENV variables for Petsc. Should be defined in .bashrc.
	export PETSC_DIR="$PWD/../../external/petsc"
  	export PETSC_ARCH=arch-darwin-c-debug 
  	export SLEPC_DIR="$PWD/../../external/petsc"

  	#Make Trilinos
	make -j4 trilinos | tee -a make_basic.log

	# Another cmake to find dependencies.
	cmake .. | tee -a make_basic.log

	make -j4 complete | tee -a make_basic.log
	./utopia_bench | tee -a make_basic.log
	./utopia_test | tee -a make_basic.log
	make install | tee -a make_basic.log
	make -j4 test_install | tee -a make_basic.log
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