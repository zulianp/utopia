#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic:.\n The output is logged in make_basic.log.\n"

_basic_build(){

	touch make_basic.log

	cmake .. -DUTOPIA_ENABLE_BLAS=ON -DUTOPIA_ENABLE_TRILINOS=OFF -DUTOPIA_ENABLE_PETSC=OFF -DUTOPIA_ENABLE_EXAMPLES=ON -DUTOPIA_ENABLE_TESTS=ON -DUTOPIA_ENABLE_SANITIZER=OFF -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia_basic | tee make_basic.log


	make -j$1 complete | tee -a make_basic.log
	./utopia_bench | tee -a make_basic.log
	./utopia_test | tee -a make_basic.log
	make install | tee -a make_basic.log
	make -j$1 test_install | tee -a make_basic.log
}


jobs=$@

if [[ -d build_basic ]]
then
	cd build_basic
	rm -rf *
	_basic_build $jobs
fi

if [[ ! -d build_basic ]]
then
	mkdir build_basic
	cd build_basic
	rm -rf *
	_basic_build $jobs
fi