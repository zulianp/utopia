#!/bin/bash

# All scripts inside https://bitbucket.org/zulianp/utopia/src/monotone_mg_refactor/utopia/scripts/installer need to be copied on the PWD
# Make sure to have a proper MPI installations and all the env variables defined (e.g. MPI_DIR)

export USER_DIR=$PWD
export INSTALL_DIR=$PWD/installation
export TRILINOS_DIR=$INSTALL_DIR/Trilinos/lib/cmake/Trilinos
export INTREPID2_DIR=$INSTALL_DIR/Trilinos/lib/cmake/Intrepid2
export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith
export UTOPIA_DIR=$INSTALL_DIR/utopia
export UTOPIA_FE_DIR=$INSTALL_DIR/utopia_fe

##################################################################

git clone https://github.com/trilinos/Trilinos.git \
    cd Trilinos && \
    mkdir build && \
    cd build && \
    source $USER_DIR/configure_trilinos.sh \
    make -j && make install

##################################################################

cd $USER_DIR
git clone https://zulianp@bitbucket.org/zulianp/par_moonolith.git
cd par_moonolith
git checkout development
mkdir build
cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$MOONOLITH_DIR && \
 make -j && make install

##################################################################

cd $USER_DIR
git clone https://bitbucket.org/zulianp/utopia.git
cd utopia
git checkout monotone_mg_refactor
git submodule update --init --recursive

cd utopia && mkdir build
cd build &&  \
    cmake .. -DUTOPIA_DEPENDENCIES_DIR=$INSTALL_DIR \
             -DCMAKE_INSTALL_PREFIX=$UTOPIA_DIR &&
    make petsc \
    export PETSC_DIR=$INSTALL_DIR/petsc \
    make -j complete \
    make install \


cd $USER_DIR/utopia/utopia_fe && mkdir build
cd build &&  \
   cmake .. -DUTOPIA_DEPENDENCIES_DIR=$INSTALL_DIR \
             -DCMAKE_INSTALL_PREFIX=$UTOPIA_FE_DIR \
             -DUTOPIA_ENABLE_LIBMESH=OFF \
             -DUTOPIA_ENABLE_MOONOLITH=ON \
             -DUTOPIA_ENABLE_INTREPID2=ON \
             -DUTOPIA_ENABLE_STK=ON && \
             -DIntrepid2_DIR=$INTREPID2_DIR
    make -j complete \
    make install \

