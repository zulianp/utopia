#ifdef UTOPIA_WITH_TRILINOS

#include "utopia_trilinos_Eval_Chop.hpp"
#include "utopia_Tpetra_Matrix_impl.hpp"

namespace utopia {

    void Chop<TpetraMatrix, TRILINOS>::apply(TpetraMatrix &mat, const Scalar &eps) {
        mat.transform_values(UTOPIA_LAMBDA(const Scalar &v) -> Scalar {
            return std::abs(v) < eps ? 0.0 : v; });
    }

    void ChopSmallerThan<TpetraMatrix, TRILINOS>::apply(TpetraMatrix &mat, const Scalar &eps) {
        mat.transform_values(UTOPIA_LAMBDA(const Scalar &v) -> Scalar {
            return v < eps ? 0.0 : v; });
    }

    void ChopGreaterThan<TpetraMatrix, TRILINOS>::apply(TpetraMatrix &mat, const Scalar &eps) {
        mat.transform_values(UTOPIA_LAMBDA(const Scalar &v) -> Scalar {
            return v > eps ? 0.0 : v; });
    }

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS
