#!/bin/bash
set -e
# Module setup and superlu
STAGE_DIR=$SCRATCH/code
INSTALL_DIR=$STAGE_DIR/installations
mkdir STAGE_DIR
cd $STAGE_DIR
source setup_fluya_env.sh

# Download Install SuperLU
git clone https://github.com/xiaoyeli/superlu.git
cd superlu
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/superlu
make -j6 && make install



cd $STAGE_DIR
#Utopia
# Download utopia if not present
if [[ ! -d utopia ]]
	then
		git clone --recurse-submodules https://bitbucket.org/zulianp/utopia.git
fi

cd utopia/utopia
git checkout cmake_refactor_fe
mkdir build_fluya && cd build_fluya
cmake .. -DUTOPIA_ENABLE_LOCAL_MODE=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fluya

make -j12 petsc
make -j12 trilinos
make yaml-cpp

cmake ..

cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DUTOPIA_ENABLE_LOCAL_MODE=OFF
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
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/par_moonolith -DCMAKE_BUILD_TYPE=Release -DMOONOLITH_ENABLE_BENCHMARK=OFF -DMOONOLITH_ENABLE_TESTING=OFF
make -j12
make install


export UTOPIA_DIR=$INSTALL_DIR/utopia_fluya
export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith

cd $STAGE_DIR
#UtopiaFE
cd utopia/utopia_fe
mkdir build_fluya && cd build_fluya
cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fe_fluya
make -j12
make install
