#ifndef UTOPIA_SPLIT_MATRIX_HPP
#define UTOPIA_SPLIT_MATRIX_HPP

#include "utopia_Views.hpp"

namespace utopia {

    template<class Scalar>
    UTOPIA_INLINE_FUNCTION static constexpr Scalar ramp_fun_positive(const Scalar &x) {
        return (device::abs(x) + x)/2;
    }
    template<class Scalar>
    UTOPIA_INLINE_FUNCTION static constexpr Scalar ramp_fun_negative(const Scalar &x) {
        return (device::abs(x) - x)/2;
    }

    template<typename Scalar, int Dim>
    UTOPIA_INLINE_FUNCTION static void split_matrix(
        const StaticMatrix<Scalar, Dim, Dim> &mat,
        StaticVector<Scalar, Dim> &values,
        StaticMatrix<Scalar, Dim, Dim> &vectors,
        StaticMatrix<Scalar, Dim, Dim> &negative,
        StaticMatrix<Scalar, Dim, Dim> &positive)
    {
        negative.set(0.0);
        positive.set(0.0);

        eig(mat, values, vectors);

        StaticVector<Scalar, Dim> v;

        for(int d = 0; d < Dim; ++d) {
            vectors.col(d, v);
            auto outer_v = outer(v, v);

            auto eig_p = ramp_fun_positive(values[d]);
            auto eig_n = ramp_fun_negative(values[d]);

            negative += eig_n * outer_v;
            positive += eig_p * outer_v;
        }
    }
}

#endif //UTOPIA_SPLIT_MATRIX_HPP
