#ifndef UTOPIA_PETSC_COND_HPP
#define UTOPIA_PETSC_COND_HPP
#ifdef WITH_SLEPC

#include "utopia_Cond.hpp"
#include "utopia_petsc_Slepc.hpp"

namespace utopia {

    template<class Matrix>
    class Cond<Matrix, PETSC> {
    public:
        using Scalar = UTOPIA_SCALAR(Matrix);
        using Vector = UTOPIA_VECTOR(Matrix);

        static Scalar apply(const Matrix &H)
        {
            SlepcSolver<Matrix, Vector> slepc;
            slepc.solver_type("arpack");
            slepc.portion_of_spectrum("smallest_real");
            slepc.number_of_eigenvalues(1);
            slepc.solve(H);
            
            Scalar small = 0, large = 0;

            Vector eigenvector;
            slepc.get_real_eigenpair(0, small, eigenvector);
            slepc.portion_of_spectrum("largest_real");
            slepc.get_real_eigenpair(0, large, eigenvector);

            return large/small;
        }
    };

}

#endif //WITH_SLEPC
#endif //UTOPIA_PETSC_COND_HPP
