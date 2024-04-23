#!/bin/bash
set -e

source setup_franetg_env_daint.sh

STAGE_DIR=$SCRATCH/code
INSTALL_DIR=$STAGE_DIR/installations

if [[ ! -d $STAGE_DIR ]]
	then
	mkdir $STAGE_DIR
fi

cd $STAGE_DIR


# Download Install SuperLU
if [[ ! -d superlu ]]
	then
		git clone https://github.com/xiaoyeli/superlu.git
fi
cd superlu
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/superlu
make -j10 && make install


# Define the install_dir of superlu.
export SuperLU_DIR=$INSTALL_DIR/superlu

# Download and install METIS and GKLIB needed for now.
cd $STAGE_DIR
if [[ ! -d GKlib ]]
	then
		git clone https://github.com/KarypisLab/GKlib.git
fi
cd GKlib
make config prefix=$INSTALL_DIR/gklib
make install

cd $STAGE_DIR
if [[ ! -d METIS ]]
	then
		git clone https://github.com/KarypisLab/METIS.git
fi
cd METIS
make config prefix=$INSTALL_DIR/metis gklib_path=$INSTALL_DIR/gklib
make install

export METIS_DIR=$INSTALL_DIR/metis

cd $STAGE_DIR
# Download Install adios2
if [[ ! -d ADIOS2 ]]
	then
		git clone https://github.com/ornladios/ADIOS2.git
fi
cd ADIOS2
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/adios2_gpu -DADIOS2_USE_CUDA=ON -DCMAKE_CXX_COMPILER=$TRILINOS_DIR/bin/nvcc_wrapper
make -j10 && make install

export ADIOS2_DIR=$INSTALL_DIR/adios2_gpu/lib64/cmake/adios2

cd $STAGE_DIR

# Download Install MARS with ADIOS2
# Need the nvcwrapper 

if [[ ! -d mars ]]
	then
		git clone https://bitbucket.org/zulianp/mars.git
fi

cd mars
# TEMPORARY
git checkout cmake_refactor
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/mars -DMARS_ENABLE_ADIOS2=ON -DADIOS2_DIR=$ADIOS2_DIR -DMARS_ENABLE_CUDA=ON
make -j10 && make install

export Mars_DIR=$INSTALL_DIR/mars

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
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_franetg_gpu -DUTOPIA_ENABLE_FRANETG_MODE=ON
make yaml-cpp

export YAML_CPP_DIR=$STAGE/utopia/external/yaml-cpp/lib64/cmake/yaml-cpp
cmake ..

make -j12
make install

export UTOPIA_DIR=$INSTALL_DIR/utopia_franetg_gpu

cd $STAGE_DIR
#UtopiaFE
cd utopia/utopia_fe
mkdir build_franetg && cd build_franetg
cmake .. -DUTOPIA_ENABLE_MARS=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fe_franetg_gpu -DUTOPIA_ENABLE_LIBMESH=OFF -DUTOPIA_ENABLE_STK=OFF -DMars_DIR=$INSTALL_DIR/mars -DCMAKE_CXX_COMPILER=$TRILINOS_DIR/bin/nvcc_wrapper
make -j10
make install

