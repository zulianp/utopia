# load necessary modules
# module load daint-mc
# module load nano
# module load cray-hdf5-parallel/1.8.16
# module load CMake
# module load cray-libsci
# module unload cray-petsc
# module load Boost
# module load Visit
# module load 



# # make dynamic linking the default
# export CRAYPE_LINK_TYPE=dynamic
# export CRAY_ADD_RPATH=yes


# export compiler wrapper paths
export CC=cc
export CXX=CC
export FC=ftn
export F90=ftn
export F77=ftn

# aliases
alias mpicc="cc"
alias mpicxx="CC"
alias mpifort="ftn"
alias mpif90="ftn"
alias mpif77="ftn"

# PETSc
export PKG_CONFIG_PATH=$CRAY_LIBSCI_PREFIX_DIR/lib/pkgconfig:$PKG_CONFIG_PATH


python2 ./configure \
--prefix=$PETSC_DIR \
--download-hypre=1 \
--with-ssl=0 \
--with-debugging=0 \
--with-pic=1 \
--with-shared-libraries=1 \
--with-cc=cc \
--with-cxx=CC \
--with-fc=ftn \
--download-fblaslapack=1 \
--download-metis=1 \
--download-parmetis=1 \
--download-mumps=1  \
--download-scalapack=1 \
--download-superlu_dist=1 \
  CC=cc CXX=ftn FC=ftn F77=ftn F90=ftn CFLAGS='-fPIC -fopenmp' CXXFLAGS='-fPIC -fopenmp' FFLAGS='-fPIC -fopenmp' FCFLAGS='-fPIC -fopenmp' F90FLAGS='-fPIC -fopenmp' F77FLAGS='-fPIC -fopenmp' PETSC_DIR='pwd'


export PETSC_DIR=/users/zulianp/petsc
export PETSC_ARCH=arch-linux2-c-opt
