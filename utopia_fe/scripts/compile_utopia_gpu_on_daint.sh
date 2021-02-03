module switch PrgEnv-cray PrgEnv-gnu
module load CMake/3.8.1
module load git/2.11.0
module load cray-hdf5-parallel
module load cudatoolkit
module load cuda_memtest/1.2.3
module load daint-gpu
module load ddt

rm -fr CMake*
cmake .. -DTRY_WITH_PETSC:BOOL=FALSE -DTRY_WITH_CUDA:BOOL=FALSE -DTRY_UTOPIA_WITH_TRILINOS:BOOL=TRUE -DTRY_WITH_EIGEN_3:BOOL=TRUE -DEIGEN3_INCLUDE_DIR:FILEPATH=/project/csstaff/nfadel/eigen/eigen-eigen-5a0156e40feb -DPETSC_DIR:FILEPATH=/project/csstaff/nfadel/shared_pasc/install/daint/gpu/petsc/  -DTRILINOS_DIR:FILEPATH=/project/csstaff/nfadel/shared_pasc/install/daint/gpu/trilinos/lib/cmake/Trilinos -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=RELWITHDEBINFO -DCMAKE_CXX_FLAGS='-fopenmp -Wall -std=c++11' -DCMAKE_CXX_COMPILER:FILEPATH=/users/nfadel/nvcc_wrapper -DCMAKE_INSTALL_PREFIX='/scratch/snx3000/nfadel/utopia-install-gpu' -DUTOPIA_INCLUDES='/project/csstaff/anfink/shared_pasc/install/daint/cpu/petsc/include'
