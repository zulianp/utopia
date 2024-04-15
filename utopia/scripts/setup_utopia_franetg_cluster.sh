#!/bin/bash
set -e

source setup_franetg_env_daint.sh

STAGE_DIR=$SCRATCH/code
INSTALL_DIR=$STAGE_DIR/installations
mkdir $STAGE_DIR
cd $STAGE_DIR


# Download Install SuperLU
git clone https://github.com/xiaoyeli/superlu.git
cd superlu
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/superlu
make -j6 && make install


# Define the install_dir of superlu.
export SuperLU_DIR=$INSTALL_DIR/superlu


# Download Install adios2
git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/adios2
make -j6 && make install

export ADIOS2_DIR=$INSTALL_DIR/adios2/lib/cmake/adios2


# UTOPIA
cd $STAGE_DIR
#Utopia
# Download utopia if not present
if [[ ! -d utopia ]]
	then
		git clone --recurse-submodules https://bitbucket.org/zulianp/utopia.git
fi


cd utopia/utopia
git checkout cmake_refactor_fe
mkdir build_franetg && cd build_franetg
cmake .. -DUTOPIA_INSTALL_TRILINOS=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_franetg -DUTOPIA_ENABLE_EIGER=ON -DUTOPIA_ENABLE_CUDA=ON

make -j12 trilinos
export CXX=$TRILINOS_DIR/bin/nvcc_wrapper

cd $STAGE_DIR

# Download Install MARS with ADIOS2
# Need the nvcwrapper 
git clone git clone https://bitbucket.org/zulianp/mars.git
cd mars
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/mars -DMARS_ENABLE_ADIOS2=ON -DADIOS2_DIR=$ADIOS2_DIR -DMARS_ENABLE_CUDA=ON
make -j6 && make install

export MARS_DIR=$INSTALL_DIR/mars

cd $STAGE_DIR/utopia/utopia/build_franetg

make yaml-cpp

cmake ..

cmake .. -DUTOPIA_ENABLE_TRILINOS=ON -DCMAKE_CXX_COMPILER=$TRILINOS_DIR/bin/nvccwrapper
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


export UTOPIA_DIR=$INSTALL_DIR/utopia_franetg
export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith

cd $STAGE_DIR
#UtopiaFE
cd utopia/utopia_fe
mkdir build_franetg && cd build_franetg
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fe_franetg
make -j12
make install

