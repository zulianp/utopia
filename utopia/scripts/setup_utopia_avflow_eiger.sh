#!/bin/bash
set -e
# Module setup and superlu
source setup_fluya_env_eiger.sh
# source setup_fluya_env_daint.sh

STAGE_DIR=$SCRATCH/code
INSTALL_DIR=$STAGE_DIR/installations

if [[ ! -d $STAGE_DIR ]]
	then
		mkdir $STAGE_DIR
fi

cd $STAGE_DIR
#Utopia
# Download utopia if not present
if [[ ! -d utopia ]]
	then
		git clone --recurse-submodules https://bitbucket.org/zulianp/utopia.git
fi

cd utopia/utopia
git checkout cmake_refactor_fe
if [[ ! -d build_avflow ]]
	then
		mkdir build_avflow
fi
cd build_avflow
cmake .. -DUTOPIA_ENABLE_LOCAL_MODE=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_avflow -DUTOPIA_ENABLE_CLUSTER=ON

make -j12 petsc
make yaml-cpp

cmake ..

cmake .. -DUTOPIA_ENABLE_AVFLOW_MODE=ON -DUTOPIA_ENABLE_LOCAL_MODE=OFF
make -j12
make install

# ParMoonolith
cd $STAGE_DIR
# Download ParMoonolith if not present
if [[ ! -d par_moonolith ]]
	then
		git clone https://bitbucket.org/zulianp/par_moonolith.git
fi
cd par_moonolith

if [[ ! -d build ]]
	then
	mkdir build
fi
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/par_moonolith -DCMAKE_BUILD_TYPE=Release -DMOONOLITH_ENABLE_BENCHMARK=OFF -DMOONOLITH_ENABLE_TESTING=OFF
make -j12
make install


export UTOPIA_DIR=$INSTALL_DIR/utopia_avflow
export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith

cd $STAGE_DIR
#UtopiaFE
cd utopia/utopia_fe

if [[ ! -d build_avflow ]]
	then
		mkdir build_avflow
fi
cd build_avflow
cmake .. -DUTOPIA_ENABLE_AVFLOW_MODE=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fe_avflow
make -j12
make install

