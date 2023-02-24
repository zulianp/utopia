#ifndef UTOPIA_TRILINOS_SOLVERS_CG_HPP
#define UTOPIA_TRILINOS_SOLVERS_CG_HPP

#include "utopia_Belos_solver.hpp"
#include "utopia_Config.hpp"
#include "utopia_ConjugateGradient.hpp"

#ifdef UTOPIA_WITH_TRILINOS_BELOS

namespace utopia {
    /**
     * @brief     Wrapper for Trilinos CG solver provided by Belos package.
     *
     * @tparam    Matrix
     * @tparam    Vector
     */
    template <typename Matrix, typename Vector>
    class ConjugateGradient<Matrix, Vector, TRILINOS> final : public BelosSolver<Matrix, Vector> {
        using Super = BelosSolver<Matrix, Vector>;
        using ParamKey = typename Super::Param::Key;

    public:
        ConjugateGradient() : Super(solver_type, solver_params_) {
            // No preconditioner
        }

        ConjugateGradient(const std::string &precond, const std::string &precond_side = "right")
            : Super(solver_type, solver_params_) {
            this->set_preconditioner(precond, precond_side);
        }

        ConjugateGradient *clone() const override { return new ConjugateGradient(*this); }

    private:
        static constexpr char solver_type[] = "Block CG";
        static constexpr std::initializer_list<ParamKey> solver_params_ = {ParamKey::BLOCK_SIZE,
                                                                           ParamKey::USE_SINGLE_RED,
                                                                           ParamKey::ORTHOGONALIZATION};
    };
}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS_BELOS
#endif  // UTOPIA_TRILINOS_SOLVERS_CG_HPP