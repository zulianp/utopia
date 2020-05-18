#ifndef UTOPIA_BLAS_EVAL_KRONECKERPRODUCT_HPP
#define UTOPIA_BLAS_EVAL_KRONECKERPRODUCT_HPP

#include "utopia_Base.hpp"
#include "utopia_Eval_KroneckerProduct.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class EvalKroneckerProduct<Matrix, Vector, BLAS> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        static void apply(const Vector &left, const Vector &right, Matrix &result) {
            const SizeType r = left.size();
            const SizeType c = right.size();

            result.resize(r, c);

            for (SizeType i = 0; i < r; ++i) {
                const auto x = left.get(i);
                for (SizeType j = 0; j < c; ++j) {
                    const auto y = right.get(j);
                    result.set(i, j, x * y);
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_BLAS_EVAL_KRONECKERPRODUCT_HPP
