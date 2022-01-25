#ifndef UTOPIA_OMNI_MATRIX_FREE_LINEAR_SOLVER_HPP
#define UTOPIA_OMNI_MATRIX_FREE_LINEAR_SOLVER_HPP

#include <memory>
#include "utopia_Describable.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Vector>
    class OmniMatrixFreeLinearSolver : public MatrixFreeLinearSolver<Vector>, public Describable {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Super = utopia::MatrixFreeLinearSolver<Vector>;
        using Super::update;

        inline void describe(std::ostream &os) const override { os << "OmniMatrixFreeLinearSolver"; }

    public:
        OmniMatrixFreeLinearSolver();
        ~OmniMatrixFreeLinearSolver() override;
        OmniMatrixFreeLinearSolver *clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void read(Input &in) override;
        void update(const Operator<Vector> &op) override;

        bool solve(const Operator<Vector> &op, const Vector &rhs, Vector &sol) override;

        void init_solver(const std::string &method, const std::vector<std::string> status_variables) override;

        void exit_solver(const SizeType &it, const Scalar &convergence_reason) override;

        bool check_convergence(const SizeType &it,
                               const Scalar &norm_grad,
                               const Scalar &rel_norm_grad,
                               const Scalar &norm_step) override;

    private:
        std::unique_ptr<MatrixFreeLinearSolver<Vector>> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_OMNI_MATRIX_FREE_LINEAR_SOLVER_HPP