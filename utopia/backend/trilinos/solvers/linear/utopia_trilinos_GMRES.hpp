#ifndef UTOPIA_TRILINOS_GMRES_HPP
#define UTOPIA_TRILINOS_GMRES_HPP

#include "utopia_Belos_solver.hpp"
#include "utopia_Config.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia {
    /**
     * @brief     Wrapper for Trilinos GMRES solver provided by Belos package.
     *
     * @tparam    Matrix
     * @tparam    Vector
     */
    template <typename Matrix, typename Vector>
    class GMRES<Matrix, Vector, TRILINOS> final : public BelosSolver<Matrix, Vector> {
        using Super = BelosSolver<Matrix, Vector>;
        using ParamKey = typename Super::Param::Key;

    public:
        GMRES() : Super(solver_type, solver_params_) {
            // No preconditioner
        }

        GMRES(const std::string &precond, const std::string &precond_side = "right")
            : Super(solver_type, solver_params_) {
            this->set_preconditioner(precond, precond_side);
        }

        GMRES *clone() const override { return new GMRES(*this); }

    private:
        static constexpr char solver_type[] = "Block GMRES";
        static constexpr std::initializer_list<ParamKey> solver_params_ = {ParamKey::BLOCK_SIZE,
                                                                           ParamKey::MAX_RESTARTS,
                                                                           ParamKey::NUM_BLOCKS,
                                                                           ParamKey::ORTHOGONALIZATION};
    };
}  // namespace utopia

#endif  // UTOPIA_TRILINOS_GMRES_HPP
