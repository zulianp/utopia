#ifndef UTOPIA_SOLVER_LAPACK_H
#define UTOPIA_SOLVER_LAPACK_H

#include "utopia_blas_Types.hpp"

#include "utopia_DirectSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia {

    namespace internals {
        bool lapack_dgesv_solve(const BlasMatrixd &A, const BlasVectord&b, BlasVectord&x);
    }

    template<>
    class LUDecomposition<BlasMatrixd, BlasVectord, BLAS> : public DirectSolver<BlasMatrixd, BlasVectord> {
    public:
        inline bool apply(const BlasVectord &b, BlasVectord &x) override
        {
            return internals::lapack_dgesv_solve(*this->get_operator(), b, x);
        }

        LUDecomposition * clone() const override
        {
            return new LUDecomposition();
        }

    };
}

#endif //UTOPIA_SOLVER_LAPACK_H
