

#ifndef UTOPIA_UTOPIA_BLASTRAITS_HPP
#define UTOPIA_UTOPIA_BLASTRAITS_HPP

#include "utopia_Traits.hpp"
#include "utopia_CRSMatrix.hpp"
#include "utopia_CCSMatrix.hpp"

namespace utopia {
    template<typename T>
    class BLASTraits {
    public:
        typedef T Scalar;
        typedef utopia::Matrix<T> Matrix;
        typedef std::vector<T> Vector;
        typedef utopia::CRSMatrix<T> CRSMatrix;
        //default sparse matrix
        typedef utopia::CRSMatrix<T> SparseMatrix;

        typedef utopia::CCSMatrix<T> CCSMatrix;
        typedef typename std::vector<T>::size_type SizeType;
        
        enum {
            Backend = BLAS
        };
    };

    UTOPIA_MAKE_TRAITS_TPL_1(std::vector, BLASTraits);
    UTOPIA_MAKE_TRAITS_DENSE_TPL_1(Matrix, BLASTraits);
    UTOPIA_MAKE_TRAITS_SPARSE_TPL_1(CRSMatrix, BLASTraits);
    UTOPIA_MAKE_TRAITS_SPARSE_TPL_1(CCSMatrix, BLASTraits);
}

#endif //UTOPIA_UTOPIA_BLASTRAITS_HPP