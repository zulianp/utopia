#ifndef UTOPIA_LINEAR_SOLVERS_INTERFACES_HPP
#define UTOPIA_LINEAR_SOLVERS_INTERFACES_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    // template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    // class BiCGStab {
    // public:
    // 	BiCGStab() {
    // 		static_assert(Backend < HOMEMADE, "BiCGStab not implemented for this backend");
    // 	}
    // };

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class GMRES {
    public:
        GMRES() { static_assert(Backend < HOMEMADE, "GMRES not implemented for this backend"); }
    };

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class MINRES {
    public:
        MINRES() { static_assert(Backend < HOMEMADE, "MINRES not implemented for this backend"); }
    };

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class SOR {
    public:
        SOR() { static_assert(Backend < HOMEMADE, "MINRES not implemented for this backend"); }
    };

    template <typename Matrix, typename Vector, int Backend = Traits<Vector>::Backend>
    class Factorization {
    public:
        Factorization() { static_assert(Backend < HOMEMADE, "Factorization not implemented for this backend"); }
    };

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class LUDecomposition {
    public:
        LUDecomposition() { static_assert(Backend < HOMEMADE, "LUDecomposition not implemented for this backend"); }
    };

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class CholeskyDecomposition {
    public:
        CholeskyDecomposition() {
            static_assert(Backend < HOMEMADE, "CholeskyDecomposition not implemented for this backend");
        }
    };

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class RedundantLinearSolver {
    public:
        RedundantLinearSolver() {
            static_assert(Backend < HOMEMADE, "RedundantLinearSolver not implemented for this backend");
        }
    };

}  // namespace utopia

#endif  // UTOPIA_LINEAR_SOLVERS_INTERFACES_HPP
