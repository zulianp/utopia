#ifndef UTOPIA_CUDA_ERROR_H
#define UTOPIA_CUDA_ERROR_H

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

namespace utopia {
    class CUDAError {
    public:
        static bool Check(const cudaError_t code);
        static bool CheckLastError();
    };
}

#endif //UTOPIA_CUDA_ERROR_H
