#ifndef UTOPIA_NORMALIZE_ROWS_HPP
#define UTOPIA_NORMALIZE_ROWS_HPP

#include "utopia_Base.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class NormalizeRows {
    public:
        using Vector = typename Traits<Matrix>::Vector;
        using Scalar = UTOPIA_SCALAR(Matrix);
        using SizeType = UTOPIA_SIZE_TYPE(Matrix);

        static void apply(Matrix &in_out, const Scalar tol) {
            Vector row_sum = sum(in_out, 1);

            each_transform(row_sum, row_sum, [tol](const SizeType i, const Scalar val) -> Scalar {
                UTOPIA_UNUSED(i);

                if (std::abs(val) > tol) {
                    return 1. / val;
                } else {
                    return 0.0;
                }
            });

            Matrix s = diag(row_sum);
            in_out = s * in_out;
        }
    };

    template <class Matrix>
    void normalize_rows(Matrix &in_out, const UTOPIA_SCALAR(Matrix) tol = 1e-15) {
        NormalizeRows<Matrix>::apply(in_out, tol);
    }

}  // namespace utopia

#endif  // UTOPIA_NORMALIZE_ROWS_HPP
