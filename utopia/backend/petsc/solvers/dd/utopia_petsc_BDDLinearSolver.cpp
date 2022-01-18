#include "utopia_petsc_BDDLinearSolver.hpp"
#include "utopia_petsc_BDDOperator.hpp"

#include "utopia_ConjugateGradient.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_OmniMatrixFreeLinearSolver_impl.hpp"
#include "utopia_petsc_LinearSolverFactory.hpp"

#include "utopia_InputParameters.hpp"

namespace utopia {

    template <typename Matrix, typename Vector>
    class BDDLinearSolver<Matrix, Vector>::Impl : public Configurable {
    public:
        void read(Input &in) override {
            in.get("use_preconditioner", use_preconditioner);
            op.read(in);

            if (solver) {
                in.get("inner_solver", *solver);
            }
        }

        Impl()
            : solver(
                  // std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>()
                  std::make_shared<OmniMatrixFreeLinearSolver<Vector>>()) {}

        std::shared_ptr<MatrixFreeLinearSolver> solver;
        BDDOperator<Matrix, Vector> op;
        bool use_preconditioner{true};
    };

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::verbose(const bool &val) {
        Super::verbose(val);
        // if (impl_->solver) {
        //     impl_->solver->verbose(val);
        // }
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::atol(const Scalar &val) {
        Super::atol(val);
        // if (impl_->solver) {
        //     impl_->solver->atol(val);
        // }
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::rtol(const Scalar &val) {
        Super::rtol(val);
        // if (impl_->solver) {
        //     impl_->solver->rtol(val);
        // }
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::stol(const Scalar &val) {
        Super::stol(val);
        // if (impl_->solver) {
        //     impl_->solver->stol(val);
        // }
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::read(Input &in) {
        Super::read(in);
        impl_->read(in);

        // if (this->verbose() && impl_->solver) {
        //     InputParameters params;
        //     params.set("verbose", true);
        //     impl_->solver->read(params);
        // }
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::set_inner_solver(const std::shared_ptr<MatrixFreeLinearSolver> &solver) {
        impl_->solver = solver;
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::init_memory(const Layout &layout) {
        assert(impl_->solver);

        if (layout.comm().size() == 1) {
            impl_->solver->init_memory(layout);
        } else {
            impl_->solver->init_memory(impl_->op.vector_layout());
        }
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::print_usage(std::ostream &os) const {
        os << "TODO\n";
    }

    template <typename Matrix, typename Vector>
    bool BDDLinearSolver<Matrix, Vector>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BDDLinearSolver::apply");

        if (this->verbose()) {
            b.comm().root_print("BDDLinearSolver::apply");
        }

        bool ok = true;
        if (b.comm().size() == 1) {
            ok = impl_->solver->solve(*this->get_operator(), b, x);
        } else {
            impl_->op.initialize(make_ref(b));

            Vector x_G;
            impl_->op.create_vector(x_G);
            ok = impl_->solver->solve(impl_->op, impl_->op.righthand_side(), x_G);

            impl_->op.finalize(x_G, x);
        }

        UTOPIA_TRACE_REGION_END("BDDLinearSolver::apply");
        return ok;
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::update(const std::shared_ptr<const Matrix> &A) {
        UTOPIA_TRACE_REGION_BEGIN("BDDLinearSolver::update");

        Super::update(A);

        if (A->comm().size() == 1) {
            // TODO
        } else {
            impl_->op.initialize(A);

            if (impl_->solver && impl_->use_preconditioner) {
                impl_->solver->set_preconditioner(impl_->op.create_preconditioner());
            }
        }

        init_memory(row_layout(*A));

        UTOPIA_TRACE_REGION_END("BDDLinearSolver::update");
    }

    template <typename Matrix, typename Vector>
    BDDLinearSolver<Matrix, Vector> *BDDLinearSolver<Matrix, Vector>::clone() const {
        auto c = utopia::make_unique<BDDLinearSolver>();
        c->copy(*this);
        return c.release();
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::copy(const BDDLinearSolver &other) {
        impl_->solver = std::shared_ptr<MatrixFreeLinearSolver>(other.impl_->solver->clone());
    }

    template <typename Matrix, typename Vector>
    BDDLinearSolver<Matrix, Vector>::BDDLinearSolver(const BDDLinearSolver &other) : Super(other) {
        copy(other);
    }

    template <typename Matrix, typename Vector>
    BDDLinearSolver<Matrix, Vector>::~BDDLinearSolver() = default;

    template <typename Matrix, typename Vector>
    BDDLinearSolver<Matrix, Vector>::BDDLinearSolver() : impl_(utopia::make_unique<Impl>()) {}

    template class BDDLinearSolver<PetscMatrix, PetscVector>;

}  // namespace utopia
