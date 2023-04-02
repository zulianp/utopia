#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script:\n"


# Should check env for trilinos_dir but for now leave like this.
_fluya_mode(){
	cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DTrilinos_DIR=/Users/dylan/Documents/Summer-Internship/Installations/Trilinos/lib/cmake/Trilinos
	make -j
	make install
	
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