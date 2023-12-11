#!/bin/bash
set -e
set -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script franetg mode:\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2


# Should check env for trilinos_dir but for now leave like this.
_franetg_mode(){
	cmake .. -DUTOPIA_ENABLE_MARS=ON -DCMAKE_INSTALL_PREFIX=$prefix

	read -p "Continue ? y/n" -n 1 -r
	echo    # (optional) move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		make -j$N_THREADS | tee make_all.log
		make install | tee make_all.log
		make -j$N_THREADS test_install | tee make_all.log
	fi
}
if [[ -d build_franetg ]]
then
	read -p "Build folder already exists. Do you want to delete and start over? y/n" -n 1 -r
	echo 
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		cd build_franetg
		rm -rf *
		_franetg_mode
	else
		cd build_franetg
		_franetg_mode
	fi
else
	if [[ ! -d build_franetg ]]
	then
		mkdir build_franetg
		cd build_franetg
		_franetg_mode
	fi
fi