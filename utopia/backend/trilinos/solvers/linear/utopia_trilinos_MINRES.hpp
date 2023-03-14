#ifndef UTOPIA_TRILINOS_MINRES_HPP
#define UTOPIA_TRILINOS_MINRES_HPP

#include "utopia_Belos_solver.hpp"
#include "utopia_Config.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia {
    /**
     * @brief     Wrapper for Trilinos MINRES solver provided by Belos package.
     *
     * @tparam    Matrix
     * @tparam    Vector
     */
    template <typename Matrix, typename Vector>
    class MINRES<Matrix, Vector, TRILINOS> final : public BelosSolver<Matrix, Vector> {
        using Super = BelosSolver<Matrix, Vector>;
        using ParamKey = typename Super::Param::Key;

    public:
        MINRES() : Super(solver_type, solver_params_) {
            // No preconditioner
        }

        MINRES(const std::string &precond, const std::string &precond_side = "right")
            : Super(solver_type, solver_params_) {
            this->set_preconditioner(precond, precond_side);
        }

        MINRES *clone() const override { return new MINRES(*this); }

    private:
        static constexpr char solver_type[] = "MINRES";
        static constexpr std::initializer_list<ParamKey> solver_params_ = {ParamKey::BLOCK_SIZE,
                                                                           ParamKey::PC_SIDE,
                                                                           ParamKey::PC_TYPE};
    };
}  // namespace utopia

#endif  // UTOPIA_TRILINOS_MINRES_HPP
