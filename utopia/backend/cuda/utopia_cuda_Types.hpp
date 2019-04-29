#ifndef UTOPIA_CUDA_TYPES_HPP
#define UTOPIA_CUDA_TYPES_HPP

#include "utopia_cuda_Traits.hpp"
#include "utopia_cuda_Matrix.hpp"
#include "utopia_cuda_Vector.hpp"
#include "utopia_cuda_Backend.hpp"
#include "utopia_Wrapper.hpp"

namespace utopia {
        ///////////////////////////////////////////////////////////////////////////////
    class CUMatrixd : public Wrapper<CUDAMatrix<double>, 2> {
    public:
        typedef utopia::Wrapper<CUDAMatrix<double>, 2> super;
        using super::super;
        using super::operator=;
    };


    class CUVectord : public Wrapper<CUDAVector<double>, 1> {
    public:
        typedef utopia::Wrapper<CUDAVector<double>, 1> super;
        using super::super;
        using super::operator=;
    };

        ///////////////////////////////////////////////////////////////////////////////

    inline void disp(const Wrapper<CUDAMatrix<double>, 2> &w) {
        backend(w).describe(w.implementation());
    }

    inline void disp(const Wrapper<CUDAVector<double>, 1> &w) {
        backend(w).describe(w.implementation());
    }

        ///////////////////////////////////////////////////////////////////////////////

    inline CUDAMatrix<double> &raw_type(Wrapper<CUDAMatrix<double>, 2> &utopiaType)
    {
        return utopiaType.implementation();
    }

    inline CUDAVector<double> &raw_type(Wrapper<CUDAVector<double>, 1> &utopiaType)
    {
        return utopiaType.implementation();
    }

        ///////////////////////////////////////////////////////////////////////////////

    inline const CUDAMatrix<double> &raw_type(const Wrapper<CUDAMatrix<double>, 2> &utopiaType)
    {
        return utopiaType.implementation();
    }

    inline const CUDAVector<double> &raw_type(const Wrapper<CUDAVector<double>, 1> &utopiaType)
    {
        return utopiaType.implementation();
    }
}

#endif //UTOPIA_CUDA_TYPES_HPP

