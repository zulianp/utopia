#ifndef UTOPIA_TRILINOS_BICGSTAB_HPP
#define UTOPIA_TRILINOS_BICGSTAB_HPP

#include "utopia_Belos_solver.hpp"
#include "utopia_BiCGStab.hpp"
#include "utopia_Config.hpp"

namespace utopia {
    /**
     * @brief     Wrapper for Trilinos BiCGStab solver provided by Belos package.
     *
     * @tparam    Matrix
     * @tparam    Vector
     */
    template <typename Matrix, typename Vector>
    class BiCGStab<Matrix, Vector, TRILINOS> final : public BelosSolver<Matrix, Vector> {
        using Super = BelosSolver<Matrix, Vector>;
        using ParamKey = typename Super::Param::Key;

    public:
        BiCGStab() : Super(solver_type) {
            // No preconditioner
        }

        BiCGStab(const std::string &precond, const std::string &precond_side = "right") : Super(solver_type) {
            this->set_preconditioner(precond, precond_side);
        }

        BiCGStab *clone() const override { return new BiCGStab(*this); }

    private:
        static constexpr char solver_type[] = "BiCGStab";
    };
}  // namespace utopia

#endif  // UTOPIA_TRILINOS_BICGSTAB_HPP
