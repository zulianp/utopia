#!/bin/bash
set -e
set -o pipefail

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Fluya mode:\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2


# Should check env for trilinos_dir but for now leave like this.
_fluya_mode(){
	touch make_fluya.log

	cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=OFF -DCMAKE_INSTALL_PREFIX=$prefix | tee make_fluya.log

	read -p "Continue ? y/n" -n 1 -r
	echo    # (optional) move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]]
	then

		make -j$N_THREADS | tee -a make_fluya.log
		make $1 install | tee -a make_fluya.log
		./utopia_bench | tee -a make_fluya.log
		./utopia_test | tee -a make_fluya.log
		make -j$N_THREADS test_install | tee -a make_fluya.log

	fi
}

if [[ -d build_fluya ]]
then
	cd build_fluya
	rm -rf *
	_fluya_mode
else
	if [[ ! -d build_fluya ]]
	then
		mkdir build_fluya
		cd build_fluya
		rm -rf *
		_fluya_mode
	fi
fi
