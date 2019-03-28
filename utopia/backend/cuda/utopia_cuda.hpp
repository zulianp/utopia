#ifndef UTOPIA_CUDA_HPP
#define UTOPIA_CUDA_HPP

#include <thrust/device_vector.h>

void run_stupid();
void run_stupid_2(thrust::device_vector<double> &vec);

#include "utopia_cuda_Backend.hpp"
#include "utopia_cuda_Matrix.hpp"
#include "utopia_cuda_Traits.hpp"
#include "utopia_cuda_Vector.hpp"
#include "utopia_cuda_Types.hpp"

#endif //UTOPIA_CUDA_HPP