#ifndef UTOPIA_PETSC_EVAL_COND_IMPL_HPP
#define UTOPIA_PETSC_EVAL_COND_IMPL_HPP

#include "utopia_Base.hpp"

#ifdef WITH_SLEPC

#include "utopia_petsc_Eval_Cond.hpp"
#include "utopia_petsc_Slepc.hpp"

namespace utopia {

    template<class Matrix>
    UTOPIA_SCALAR(Matrix) Cond<Matrix, PETSC>::apply(const Matrix &H)
    {
        using Vector = utopia::Wrapper<UTOPIA_VECTOR(Matrix), 1>;

        Scalar small = 0, large = 0;

        Vector eigenvector;
        
        SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> slepc;
        // slepc.solver_type("lapack");
        // slepc.solver_type("arpack");
        slepc.number_of_eigenvalues(1);

        //small eig
        slepc.portion_of_spectrum("smallest_real");
        slepc.solve(H);        
        slepc.get_real_eigenpair(0, small, eigenvector);
        
        //large eig
        slepc.portion_of_spectrum("largest_real");
        slepc.solve(H);
        slepc.get_real_eigenpair(0, large, eigenvector);

        // std::cout << "small: " << small << ", large: " << large << std::endl;
        return large/small;
    }
}

#endif //WITH_SLEPC
#endif //UTOPIA_PETSC_EVAL_COND_IMPL_HPP
