#ifndef UTOPIA_PETSC_EVAL_COND_IMPL_HPP
#define UTOPIA_PETSC_EVAL_COND_IMPL_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_SLEPC

#include "utopia_petsc_Eval_Cond.hpp"
#include "utopia_petsc_Slepc.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    template <class Matrix>
    UTOPIA_SCALAR(Matrix)
    Cond<Matrix, PETSC>::apply(const Matrix &H) {
        using Vector = typename Traits<Matrix>::Vector;

        Scalar small = 0, large = 0;

        Vector eigenvector;

        SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> slepc;

#ifdef SLEPC_HAVE_ARPACK
        slepc.solver_type("arpack");
#else
#ifdef SLEPC_HAVE_PRIMME
        slepc.solver_type("primme");
#else
#ifdef SLEPC_HAVE_TRLAN
        slepc.solver_type("trlan");
#endif
#endif
#endif

        slepc.problem_type("non_hermitian");
        slepc.number_of_eigenvalues(1);

        // small eig
        slepc.portion_of_spectrum("smallest_real");
        slepc.max_it(40000);
        // slepc.portion_of_spectrum("smallest_magnitude");
        slepc.solve(H);
        slepc.get_real_eigenpair(0, small, eigenvector);

        // large eig
        slepc.portion_of_spectrum("largest_real");
        // slepc.portion_of_spectrum("largest_magnitude");
        slepc.solve(H);
        slepc.get_real_eigenpair(0, large, eigenvector);

        const Scalar result = large / small;
        // utopia::out() <<"small: " << small << ", large: " << large << ", cond: " << result << std::endl;
        return result;
    }
}  // namespace utopia

#endif  // UTOPIA_ENABLE_SLEPC
#endif  // UTOPIA_PETSC_EVAL_COND_IMPL_HPP
