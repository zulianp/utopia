#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Franetg mode:\n"


# Should check env for trilinos_dir but for now leave like this.
_fluya_mode(){
	cmake .. -DUTOPIA_ENABLE_TRILINOS=ON -DUTOPIA_ENABLE_MARS=ON -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia
	make -j4
	./utopia_test
	./utopia_bench
	make install
	make test_install
	
}

if [[ -d build_fluya ]]
then
	cd build_fluya
	rm -rf *
	_fluya_mode
fi

if [[ ! -d build_fluya ]]
then
	mkdir build_fluya
	cd build_fluya
	rm -rf *
	_fluya_mode
fi