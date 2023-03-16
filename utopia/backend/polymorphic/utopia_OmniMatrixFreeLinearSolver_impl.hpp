#ifndef UTOPIA_OMNI_MATRIX_FREE_LINEAR_SOLVER_IMPL_HPP
#define UTOPIA_OMNI_MATRIX_FREE_LINEAR_SOLVER_IMPL_HPP

#include "utopia_MatrixFreeLinearSolverFactory.hpp"
#include "utopia_OmniMatrixFreeLinearSolver.hpp"

namespace utopia {

    template <class Vector>
    void OmniMatrixFreeLinearSolver<Vector>::read(Input &in) {
        Super::read(in);

        std::string backend = Traits<Vector>::backend_info().get_name();
        std::string type = "gmres";

        in.get("backend", backend);
        in.get("type", type);

        impl_ = MatrixFreeLinearSolverFactory<Vector>::new_linear_solver(type);
        impl_->read(in);
    }

    template <class Vector>
    OmniMatrixFreeLinearSolver<Vector>::OmniMatrixFreeLinearSolver()
        : impl_(MatrixFreeLinearSolverFactory<Vector>::default_linear_solver()) {}

    template <class Vector>
    OmniMatrixFreeLinearSolver<Vector>::~OmniMatrixFreeLinearSolver() = default;

    template <class Vector>
    OmniMatrixFreeLinearSolver<Vector> *OmniMatrixFreeLinearSolver<Vector>::clone() const {
        auto cloned = utopia::make_unique<OmniMatrixFreeLinearSolver>();

        if (this->impl_) {
            cloned->impl_ = std::unique_ptr<MatrixFreeLinearSolver<Vector>>(this->impl_->clone());
        }

        return cloned.release();
    }

    template <class Vector>
    bool OmniMatrixFreeLinearSolver<Vector>::apply(const Vector &rhs, Vector &sol) {
        assert(static_cast<bool>(impl_));
        if (!impl_) return false;

        return impl_->apply(rhs, sol);
    }

    template <class Vector>
    void OmniMatrixFreeLinearSolver<Vector>::update(const Operator<Vector> &op) {
        assert(static_cast<bool>(impl_));

        // FIXME
        // Super::update(op);

        impl_->update(op);
    }

    template <class Vector>
    bool OmniMatrixFreeLinearSolver<Vector>::solve(const Operator<Vector> &op, const Vector &rhs, Vector &sol) {
        assert(static_cast<bool>(impl_));
        return impl_->solve(op, rhs, sol);
    }

    template <class Vector>
    void OmniMatrixFreeLinearSolver<Vector>::init_solver(const std::string &method,
                                                         const std::vector<std::string> status_variables) {
        assert(static_cast<bool>(impl_));
        impl_->init_solver(method, status_variables);
    }

    template <class Vector>
    void OmniMatrixFreeLinearSolver<Vector>::exit_solver(const SizeType &it, const Scalar &convergence_reason) {
        assert(static_cast<bool>(impl_));
        impl_->exit_solver(it, convergence_reason);
    }

    template <class Vector>
    bool OmniMatrixFreeLinearSolver<Vector>::check_convergence(const SizeType &it,
                                                               const Scalar &norm_grad,
                                                               const Scalar &rel_norm_grad,
                                                               const Scalar &norm_step) {
        assert(static_cast<bool>(impl_));
        return impl_->check_convergence(it, norm_grad, rel_norm_grad, norm_step);
    }

}  // namespace utopia

#endif  // UTOPIA_OMNI_MATRIX_FREE_LINEAR_SOLVER_IMPL_HPP
