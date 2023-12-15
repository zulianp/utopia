#ifndef UTOPIA_TRILINOS_FACTORIZATION_HPP
#define UTOPIA_TRILINOS_FACTORIZATION_HPP

#include "utopia_Amesos2_solver.hpp"
#include "utopia_Config.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia {
#ifdef UTOPIA_WITH_TRILINOS_AMESOS2
    /**
     * @brief     Wrapper for Trilinos Factorization solver provided by Amesos2 package.
     *
     * @tparam    Matrix
     * @tparam    Vector
     */
    template <typename Matrix, typename Vector>
    class Factorization<Matrix, Vector, TRILINOS> final : public Amesos2Solver<Matrix, Vector> {
        using Super = Amesos2Solver<Matrix, Vector>;

    public:
        Factorization() : Super("KLU2", {{"trans", {"Trans", "One of NOTRANS, TRANS, CONJ", "NOTRANS"}}}) {}

        auto clone() const -> Factorization* override { return new Factorization(*this); }
    };
#endif  // UTOPIA_WITH_TRILINOS_AMESOS2
}  // namespace utopia

#endif  // UTOPIA_TRILINOS_FACTORIZATION_HPP
