#include "utopia_cuda_Error.hpp"

namespace utopia {
    bool CUDAError::Check(const cudaError_t code)
    {
        if (code != cudaSuccess)
        {
            fprintf(stderr,"CUDA Error: %s\n", cudaGetErrorString(code));
            return false;
        }

        return true;
    }

    bool CUDAError::CheckLastError()
    {
        return Check(cudaPeekAtLastError());
    }
}
