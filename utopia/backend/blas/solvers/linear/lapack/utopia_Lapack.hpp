//
// Created by Alessandro Rigazzi on 22/05/15.
//

#ifndef UTOPIA_SOLVER_LAPACK_H
#define UTOPIA_SOLVER_LAPACK_H

#include "utopia_blas_Types.hpp"

#include "utopia_DirectSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia {

    namespace internals {
        bool lapack_dgesv_solve(const Matrixd &A, const Vectord&b, Vectord&x);
    }

    template<>
    class LUDecomposition<Matrixd, Vectord, BLAS> : public DirectSolver<Matrixd, Vectord> {
    public:
        inline bool apply(const Vectord &b, Vectord &x) override
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
