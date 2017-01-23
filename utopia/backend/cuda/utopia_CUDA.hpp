#ifndef UTOPIA_CUDA_HPP
#define UTOPIA_CUDA_HPP 

#include <thrust/device_vector.h>

void run_stupid();
void run_stupid_2(thrust::device_vector<double> &vec);

#include "utopia_CUDABackend.hpp"
#include "utopia_CUDAMatrix.hpp"
#include "utopia_CUDATraits.hpp"
#include "utopia_CUDAVector.hpp"
#include "utopia_CUDATypes.hpp"

#endif //UTOPIA_CUDA_HPP