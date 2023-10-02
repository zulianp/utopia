#!/bin/bash

today=$(date)
printf "%s\n" "$today"
printf "Testing Cmake Script Fluya mode:\n"


# Should check env for trilinos_dir but for now leave like this.
_fluya_mode(){
	touch make_fluya.log

	cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DUTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL=OFF -DCMAKE_INSTALL_PREFIX=/Users/dylan/Documents/Summer-Internship/Installations/utopia | tee make_fluya.log

	make -j6 | tee -a make_fluya.log
	make install | tee -a make_fluya.log
	./utopia_bench | tee -a make_fluya.log
	./utopia_test | tee -a make_fluya.log
	make -j4 test_install | tee -a make_fluya.log
}

# if [[ -d build_fluya ]]
# then
# 	cd build_fluya
# 	rm -rf *
# 	_fluya_mode
# fi

if [[ ! -d build_fluya ]]
then
	mkdir build_fluya
	cd build_fluya
	rm -rf *
	_fluya_mode
fi