#ifndef UTOPIA_BACKENDS_HPP
#define UTOPIA_BACKENDS_HPP 

#include "utopia_Base.hpp"

#ifdef WITH_BLAS
#include "utopia_BLAS.hpp"
#endif //WITH_BLAS

#ifdef WITH_PETSC
#include "utopia_PETSc.hpp"
#endif //WITH_PETSC

#ifdef WITH_CUDA
#include "utopia_CUDA.hpp"
#endif //WITH_CUDA

#ifdef WITH_UTOPIA_OPENCL
#include "utopia_OpenCL.hpp"
#endif //WITH_UTOPIA_OPENCL

#endif //UTOPIA_BACKENDS_HPP
