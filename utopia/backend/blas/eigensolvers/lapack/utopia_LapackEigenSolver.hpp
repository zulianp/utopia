#ifndef UTOPIA_LAPACK_EIGEN_SOLVER_HPP
#define UTOPIA_LAPACK_EIGEN_SOLVER_HPP

#include "utopia_blas_Types.hpp"

namespace utopia {
    class LapackEigenSolver {
    public:
        ///@brief Solves A * x = lambda * B * x
        bool spd_eig(const BlasMatrixd &A, BlasVectord &evalues, BlasMatrixd &evectors) const;
        bool spd_geig(const BlasMatrixd &A, const BlasMatrixd &B, BlasVectord &evalues, BlasMatrixd &evectors) const;
        bool spd_geig_small(const BlasMatrixd &A, const BlasMatrixd &B, const double upper_bound, BlasVectord &evalues, BlasMatrixd &evectors) const;

    private:
        bool fix_sizes(const SizeType m, BlasVectord &evalues, BlasMatrixd &evectors) const;
    };

    inline bool spd_geig(const BlasMatrixd &A, const BlasMatrixd &B, BlasVectord &evalues, BlasMatrixd &evectors)
    {
        return LapackEigenSolver().spd_geig(A, B, evalues, evectors);
    }

    inline bool spd_eig(const BlasMatrixd &A, BlasVectord &evalues, BlasMatrixd &evectors)
    {
        return LapackEigenSolver().spd_eig(A, evalues, evectors);
    }

    inline bool spd_geig_small(const BlasMatrixd &A, const BlasMatrixd &B, const double upper_bound, BlasVectord &evalues, BlasMatrixd &evectors)
    {
        return LapackEigenSolver().spd_geig_small(A, B, upper_bound, evalues, evectors);
    }
}

#endif //UTOPIA_LAPACK_EIGEN_SOLVER_HPP
