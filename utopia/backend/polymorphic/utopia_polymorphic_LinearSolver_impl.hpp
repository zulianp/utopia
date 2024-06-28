#include "utopia_LinearSolverFactory.hpp"
#include "utopia_PreconditionedSolverInterface.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

namespace utopia {

#define UTOPIA_READ_ENV(name, conversion) \
    do {                                  \
        char *var = getenv(#name);        \
        if (var) {                        \
            name = conversion(var);       \
        }                                 \
    } while (0)

    template <class Matrix, class Vector>
    void OmniLinearSolver<Matrix, Vector>::read(Input &in) {
        Super::read(in);

        // std::string backend = Traits<Vector>::backend_info().get_name();
        std::string type = "gmres";

        // in.get("backend", backend);
        in.get("type", type);

        // utopia::out() << "[OmniLinearSolver] type: " << type << "\n";

        impl_ = LinearSolverFactory<Matrix, Vector>::new_linear_solver(type);

        std::string preconditioner = "";
        in.get("precond", [&](Input &node) {
            std::string preconditioner_type = "";
            node.get("type", preconditioner_type);

            // utopia::out() << "[OmniLinearSolver] preconditioner_type: " << preconditioner_type << "\n";

            auto ps = dynamic_cast<PreconditionedSolverInterface<Vector> *>(impl_.get());

            if (!ps) {
                assert(false);
                Utopia::Abort("[Error] OmniLinearSolver unable to pass preconditioner to solver!");
            }

            auto prec = LinearSolverFactory<Matrix, Vector>::new_linear_solver(preconditioner_type);
            if (prec) {
                ps->set_preconditioner(std::move(prec));
            } else {
                assert(false);
                Utopia::Abort("[Error] OmniLinearSolver (no such preconditioner)");
            }
        });

        impl_->read(in);
    }

    template <class Matrix, class Vector>
    void OmniLinearSolver<Matrix, Vector>::verbose(const bool &val) {
        Super::verbose(val);
        if (impl_) {
            impl_->verbose(val);
        }
    }

    template <class Matrix, class Vector>
    void OmniLinearSolver<Matrix, Vector>::clear() {
        if (impl_) {
            impl_->clear();
        }
    }

    template <class Matrix, class Vector>
    OmniLinearSolver<Matrix, Vector>::OmniLinearSolver()
        : impl_(LinearSolverFactory<Matrix, Vector>::default_linear_solver()) {}

    template <class Matrix, class Vector>
    OmniLinearSolver<Matrix, Vector>::~OmniLinearSolver() = default;

    template <class Matrix, class Vector>
    OmniLinearSolver<Matrix, Vector> *OmniLinearSolver<Matrix, Vector>::clone() const {
        auto cloned = utopia::make_unique<OmniLinearSolver>();

        if (this->impl_) {
            cloned->impl_ = std::unique_ptr<LinearSolver<Matrix, Vector>>(this->impl_->clone());
        }

        return cloned.release();
    }

    template <class Matrix, class Vector>
    bool OmniLinearSolver<Matrix, Vector>::apply(const Vector &rhs, Vector &sol) {
        assert(static_cast<bool>(impl_));
        if (!impl_) return false;

        return impl_->apply(rhs, sol);
    }

    template <class Matrix, class Vector>
    void OmniLinearSolver<Matrix, Vector>::update(const std::shared_ptr<const Matrix> &mat) {
        assert(static_cast<bool>(impl_));

        Super::update(mat);
        impl_->update(mat);
    }

}  // namespace utopia
