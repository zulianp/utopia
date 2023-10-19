#include "utopia_petsc_RebalancedSolver.hpp"

#ifdef UTOPIA_WITH_PARMETIS

#include "utopia_ConjugateGradient.hpp"

#include "utopia_petsc_Decompose.hpp"

#include "utopia_polymorphic_LinearSolver.hpp"

#include "utopia_petsc_Matrix_impl.hpp"

#include <vector>

namespace utopia {

    class RebalancedSolver::Impl {
    public:
        using IndexArray = Traits<PetscMatrix>::IndexArray;
        //
        std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> solver;

        PetscMatrix op;
        PetscVector rhs, sol;

        std::vector<int> partitioning, inverse_partitioning;
        IndexArray permutation, inverse_permutation;

        int block_size{1};
        bool keep_symbolic_factorization{false};

        bool empty() { return partitioning.empty(); }

        void clear() {
            op.clear();
            rhs.clear();
            sol.clear();
            partitioning.clear();
            inverse_partitioning.clear();
            permutation.clear();
            inverse_permutation.clear();
        }
    };

    void RebalancedSolver::read(Input &in) {
        Super::read(in);
        in.get("inner_solver", *impl_->solver);
        in.get("block_size", impl_->block_size);
        in.get("keep_symbolic_factorization", impl_->keep_symbolic_factorization);
    }

    RebalancedSolver::RebalancedSolver() : impl_(utopia::make_unique<Impl>()) {
        set_solver(std::make_shared<OmniLinearSolver<PetscMatrix, PetscVector>>());
    }

    RebalancedSolver::~RebalancedSolver() {}

    void RebalancedSolver::set_solver(const std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> &solver) {
        impl_->solver = solver;
    }

    RebalancedSolver *RebalancedSolver::clone() const {
        auto ptr = utopia::make_unique<RebalancedSolver>();
        ptr->impl_->solver = std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>>(impl_->solver->clone());
        return ptr.release();
    }

    void RebalancedSolver::clear() {
        UTOPIA_TRACE_SCOPE("RebalancedSolver::clear");
        impl_->clear();
    }

    void RebalancedSolver::update(const std::shared_ptr<const PetscMatrix> &op) {
        UTOPIA_TRACE_SCOPE("RebalancedSolver::update");

        Super::update(op);

        if (op->comm().size() == 1) {
            impl_->solver->update(op);
        } else {
            if (impl_->empty() || !impl_->keep_symbolic_factorization) {
                if (impl_->block_size == 1) {
                    initialize_rebalance(*op,
                                         impl_->partitioning,
                                         impl_->permutation,
                                         impl_->inverse_partitioning,
                                         impl_->inverse_permutation);

                } else {
                    initialize_rebalance_block(impl_->block_size,
                                               *op,
                                               impl_->partitioning,
                                               impl_->permutation,
                                               impl_->inverse_partitioning,
                                               impl_->inverse_permutation);
                }
            }

            utopia::redistribute_from_permutation(*op, impl_->permutation, impl_->op);

            if (0) {
                std::stringstream ss;
                ss << impl_->op.local_rows();
                op->comm().synched_print(ss.str());
            }

            impl_->solver->update(make_ref(impl_->op));
        }
    }

    bool RebalancedSolver::apply(const PetscVector &rhs, PetscVector &sol) {
        UTOPIA_TRACE_SCOPE("RebalancedSolver::apply");

        if (rhs.comm().size() == 1) {
            return impl_->solver->apply(rhs, sol);
        } else {
            // redist rhs and sol
            redistribute_from_permutation(rhs, impl_->permutation, impl_->rhs);
            redistribute_from_permutation(sol, impl_->permutation, impl_->sol);

            bool ok = impl_->solver->apply(impl_->rhs, impl_->sol);

            // redistribute solution
            redistribute_from_permutation(impl_->sol, impl_->inverse_permutation, sol);
            return ok;
        }
    }

}  // namespace utopia

#endif  // UTOPIA_WITH_PARMETIS
