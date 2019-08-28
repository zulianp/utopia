#ifndef UTOPIA_PETSC_COND_HPP
#define UTOPIA_PETSC_COND_HPP

#include "utopia_Base.hpp"

#ifdef WITH_SLEPC

#include "utopia_Cond.hpp"
// #include "utopia_petsc_Slepc.hpp"

namespace utopia {

    template<class Matrix>
    class Cond<Matrix, PETSC> {
    public:
        using Scalar = UTOPIA_SCALAR(Matrix);
        // using Vector = utopia::Wrapper<UTOPIA_VECTOR(Matrix), 1>;

        static Scalar apply(const Matrix &H);
        // {
        //     SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> slepc;
        //     // slepc.solver_type("lapack");
        //     slepc.portion_of_spectrum("smallest_real");
        //     slepc.number_of_eigenvalues(1);
        //     slepc.solve(H);

        //     Scalar small = 0, large = 0;

        //     Vector eigenvector;
        //     slepc.get_real_eigenpair(0, small, eigenvector);
        //     slepc.portion_of_spectrum("largest_real");
        //     slepc.get_real_eigenpair(0, large, eigenvector);

        //     return large/small;
        // }
    };

}

#endif //WITH_SLEPC
#endif //UTOPIA_PETSC_COND_HPP
