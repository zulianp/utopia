#ifndef UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP
#define UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP


#include "utopia_Preconditioner.hpp"

namespace  utopia
{
    /**
     * @brief      The base class for matrix free linear solvers.
     * @tparam     Vector
     */
    template<class Vector>
    class MatrixFreeLinearSolver {
    public:
        virtual ~MatrixFreeLinearSolver() {}
        virtual bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) = 0;
    };
}

#endif //UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP
