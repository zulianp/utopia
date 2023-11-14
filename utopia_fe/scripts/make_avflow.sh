#!/bin/bash
set -e
set -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script avflow mode:\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2


# Should check env for trilinos_dir but for now leave like this.
_avflow_mode(){
	cmake $SCRIPT_DIR/../ -DUTOPIA_ENABLE_TRILINOS=OFF -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_LIBMESH=ON -DCMAKE_INSTALL_PREFIX=$prefix
	make -j$N_THREADS
	./utopia_test
	./utopia_bench
	make -j$N_THREADS install
	make -j$N_THREADS test_install
	
}
if [[ -d build_avflow ]]
then
	read -p "Build folder already exists. Do you want to delete and start over? y/n" -n 1 -r
	echo 
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		cd build_avflow
		rm -rf *
		_avflow_mode
	else
		cd build_avflow
		_avflow_mode
	fi
else
	if [[ ! -d build_avflow ]]
	then
		mkdir build_avflow
		cd build_avflow
		_avflow_mode
	fi
fi