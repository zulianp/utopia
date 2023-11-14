#!/bin/bash
set -e
set -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Basic:.\nThe output is logged in make_basic.log.\n"

if [[ -z "$1" ]]
then
	N_THREADS=4
else
	N_THREADS=$1
fi

prefix=$2

_basic_build(){

	touch make_basic.log

	cmake $SCRIPT_DIR/../
}

if [[ -d build_basic ]]
then
	read -p "Build folder already exists. Do you want to delete and start over? y/n" -n 1 -r
	echo 
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		cd build_basic
		rm -rf *
		_basic_build
	else
		cd build_basic
		_basic_build
	fi
else
	if [[ ! -d build_basic ]]
	then
		mkdir build_basic
		cd build_basic
		_basic_build
	fi
fi