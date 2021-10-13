#include "utopia_petsc_BDDLinearSolver.hpp"
#include "utopia_petsc_BDDOperator.hpp"

#include "utopia_ConjugateGradient.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    template <typename Matrix, typename Vector>
    class BDDLinearSolver<Matrix, Vector>::Impl : public Configurable {
    public:
        void read(Input &in) override {
            if (solver) {
                solver->read(in);
            }
        }

        Impl() : solver(std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>()) {}

        std::shared_ptr<MatrixFreeLinearSolver> solver;
        BDDOperator<Matrix, Vector> op;
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
        impl_->read(in);
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::set_inner_solver(const std::shared_ptr<MatrixFreeLinearSolver> &solver) {
        impl_->solver = solver;
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::init_memory(const Layout &layout) {
        UTOPIA_UNUSED(layout);
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::print_usage(std::ostream &os) const {
        os << "TODO\n";
    }

    template <typename Matrix, typename Vector>
    bool BDDLinearSolver<Matrix, Vector>::apply(const Vector &b, Vector &x) {
        impl_->op.initialize(make_ref(b));

        Vector x_G;
        impl_->op.create_vector(x_G);
        bool ok = impl_->solver->solve(impl_->op, impl_->op.righthand_side(), x_G);

        impl_->op.finalize(x_G, x);

        return ok;
    }

    template <typename Matrix, typename Vector>
    void BDDLinearSolver<Matrix, Vector>::update(const std::shared_ptr<const Matrix> &A) {
        impl_->op.initialize(A);
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
