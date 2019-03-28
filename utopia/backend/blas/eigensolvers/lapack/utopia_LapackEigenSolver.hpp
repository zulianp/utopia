#ifndef UTOPIA_LAPACK_EIGEN_SOLVER_HPP
#define UTOPIA_LAPACK_EIGEN_SOLVER_HPP

#include "utopia_blas_Types.hpp"

namespace utopia {
    class LapackEigenSolver {
    public:
        ///@brief Solves A * x = lambda * B * x
        bool spd_eig(const Matrixd &A, Vectord &evalues, Matrixd &evectors) const;
        bool spd_geig(const Matrixd &A, const Matrixd &B, Vectord &evalues, Matrixd &evectors) const;
        bool spd_geig_small(const Matrixd &A, const Matrixd &B, const double upper_bound, Vectord &evalues, Matrixd &evectors) const;

    private:
        bool fix_sizes(const SizeType m, Vectord &evalues, Matrixd &evectors) const;
    };

    inline bool spd_geig(const Matrixd &A, const Matrixd &B, Vectord &evalues, Matrixd &evectors)
    {
        return LapackEigenSolver().spd_geig(A, B, evalues, evectors);
    }

    inline bool spd_eig(const Matrixd &A, Vectord &evalues, Matrixd &evectors)
    {
        return LapackEigenSolver().spd_eig(A, evalues, evectors);
    }

    inline bool spd_geig_small(const Matrixd &A, const Matrixd &B, const double upper_bound, Vectord &evalues, Matrixd &evectors)
    {
        return LapackEigenSolver().spd_geig_small(A, B, upper_bound, evalues, evectors);
    }
}

#endif //UTOPIA_LAPACK_EIGEN_SOLVER_HPP
