
#ifndef UTOPIA_BLAS_TYPES_HPP
#define UTOPIA_BLAS_TYPES_HPP

#include "utopia_blas_Matrix.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_blas_CRSMatrix.hpp"

#include "utopia_Wrapper.hpp"

namespace utopia {
    typedef utopia::Wrapper<utopia::Matrix<double>, 2>    Matrixd;
    typedef utopia::Wrapper<std::vector<double>, 1>       Vectord;
    typedef utopia::Wrapper<utopia::CRSMatrix<double>, 2> CRSMatrixd;
    typedef utopia::Wrapper<utopia::CCSMatrix<double>, 2> CCSMatrixd;

    ///////////////////////////////////////////////////////////////////////////////

    inline Matrix<double> &raw_type(Wrapper<Matrix<double>, 2> &utopiaType)
    {
        return utopiaType.implementation();
    }

    inline std::vector<double> &raw_type(Wrapper<std::vector<double>, 1> &utopiaType)
    {
        return utopiaType.implementation();
    }

    ///////////////////////////////////////////////////////////////////////////////

    inline const Matrix<double> &raw_type(const Wrapper<Matrix<double>, 2> &utopiaType)
    {
        return utopiaType.implementation();
    }

    inline const std::vector<double> &raw_type(const Wrapper<std::vector<double>, 1> &utopiaType)
    {
        return utopiaType.implementation();
    }

}

#endif //UTOPIA_BLAS_TYPES_HPP
