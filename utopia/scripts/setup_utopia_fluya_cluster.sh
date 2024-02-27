#!/bin/bash
set -e
# Module setup and superlu
cd $SCRATCH
module load cray
module load PrgEnv-gnu
module load CMake
module load gcc/11.2.0  cray-mpich/8.1.12  cray-hdf5-parallel/1.12.1.3
module load cray-netcdf-hdf5parallel/4.8.1.3


# Download Install SuperLU
git clone https://github.com/xiaoyeli/superlu.git
cd superlu
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$SCRATCH/Installations/superlu
make -j6 && make install

cd $SCRATCH
#Utopia
# Download utopia if not present
if [[ ! -d utopia ]]
	then
		git clone --recurse-submodules https://bitbucket.org/zulianp/utopia.git
fi

cd utopia/utopia
git checkout cmake_refactor_fe
mkdir build_fluya && cd build_fluya
cmake .. -DUTOPIA_ENABLE_LOCAL_MODE=ON -DCMAKE_INSTALL_PREFIX=$SCRATCH/Installations/utopia_fluya -DUTOPIA_ENABLE_EIGER=ON

make -j12 petsc
make -j12 trilinos
make yaml-cpp

cmake ..

cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DUTOPIA_ENABLE_LOCAL_MODE=OFF
make -j12
make install

# ParMoonolith
cd $SCRATCH
# Download ParMoonolith if not present
if [[ ! -d par_moonolith ]]
	then
		git clone https://bitbucket.org/zulianp/par_moonolith.git
fi
cd par_moonolith
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$SCRATCH/Installations/par_moonolith -DCMAKE_BUILD_TYPE=Release -DMOONOLITH_ENABLE_BENCHMARK=OFF -DMOONOLITH_ENABLE_TESTING=OFF
make -j12
make install


export UTOPIA_DIR=$SCRATCH/Installations/utopia_fluya
export MOONOLITH_DIR=$SCRATCH/Installations/par_moonolith

cd $SCRATCH
#UtopiaFE
cd utopia/utopia_fe
mkdir build_fluya && cd build_fluya
cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DCMAKE_INSTALL_PREFIX=$SCRATCH/Installations/utopia_fe_fluya
make -j12
make install
