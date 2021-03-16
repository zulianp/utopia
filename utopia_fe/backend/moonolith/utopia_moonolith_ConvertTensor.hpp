#ifndef UTOPIA_MOONOLITH_CONVERT_TENSOR_HPP
#define UTOPIA_MOONOLITH_CONVERT_TENSOR_HPP

#include "utopia_ConvertTensor.hpp"

#include "utopia_fe_base.hpp"

namespace moonolith {

    template <typename T>
    class SparseMatrix;
}

namespace utopia {

    template <typename T>
    class Traits<::moonolith::SparseMatrix<T>> {
    public:
        static const int Order = 2;
        using Scalar = T;
    };

    template <>
    class ConvertTensor<::moonolith::SparseMatrix<UScalar>, USparseMatrix, 2> {
    public:
        static void apply(const ::moonolith::SparseMatrix<UScalar> &in, USparseMatrix &out);
    };

    template <typename T, class To, int Order>
    inline void convert(const ::moonolith::SparseMatrix<T> &from, Tensor<To, Order> &to) {
        ConvertTensor<::moonolith::SparseMatrix<T>, To, Order>::apply(from, to.derived());
    }
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_CONVERT_TENSOR_HPP