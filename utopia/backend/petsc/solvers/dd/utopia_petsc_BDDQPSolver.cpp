#include "utopia_petsc_BDDQPSolver.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_petsc_BDDQPSolver.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

#include "utopia_MPRGP.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_petsc_BDDOperator.hpp"
#include "utopia_petsc_Factorization.hpp"

#include "utopia_DeviceView.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class BDDQPSolver<Matrix, Vector, Backend>::Impl {
    public:
        BDDOperator<Matrix, Vector> op;
        std::shared_ptr<MatrixFreeQPSolver> solver;
        bool debug{false};
        bool determine_boolean_selector{false};
        Scalar infinity{1000};
        bool use_preconditioner{true};
        std::shared_ptr<Vector> boolean_selector;
    };

    template <class Matrix, class Vector, int Backend>
    BDDQPSolver<Matrix, Vector, Backend>::BDDQPSolver(std::shared_ptr<MatrixFreeQPSolver> solver)
        : impl_(utopia::make_unique<Impl>()) {
        set_solver(solver);
    }

    template <class Matrix, class Vector, int Backend>
    BDDQPSolver<Matrix, Vector, Backend>::BDDQPSolver() : impl_(utopia::make_unique<Impl>()) {
        auto mprgp = std::make_shared<MPRGP<Matrix, Vector>>();
        set_solver(mprgp);
    }

    template <class Matrix, class Vector, int Backend>
    BDDQPSolver<Matrix, Vector, Backend>::~BDDQPSolver() = default;

    template <class Matrix, class Vector, int Backend>
    BDDQPSolver<Matrix, Vector, Backend> *BDDQPSolver<Matrix, Vector, Backend>::clone() const {
        auto ptr = utopia::make_unique<BDDQPSolver>();
        if (this->has_bound()) {
            ptr->set_box_constraints(this->get_box_constraints());
        }

        return ptr.release();
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::read(Input &in) {
        QPSolver<Matrix, Vector>::read(in);
        in.get("debug", impl_->debug);
        in.get("infinity", impl_->infinity);
        in.get("determine_boolean_selector", impl_->determine_boolean_selector);
        in.get("use_preconditioner", impl_->use_preconditioner);

        if (impl_->solver) {
            // impl_->solver->read(in);
            in.get("inner_solver", *impl_->solver);
        }

        impl_->op.read(in);
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::print_usage(std::ostream &os) const {
        Super::print_usage(os);
    }

    template <class Matrix, class Vector, int Backend>
    bool BDDQPSolver<Matrix, Vector, Backend>::smooth(const Vector &, Vector &) {
        assert(false && "IMPLEMENT ME");
        return false;
    }

    template <class Matrix, class Vector, int Backend>
    bool BDDQPSolver<Matrix, Vector, Backend>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BDDQPSolver::apply");

        if (this->verbose()) {
            b.comm().root_print("BDDQPSolver::apply");
        }

        bool ok = true;
        if (b.comm().size() == 1) {
            impl_->solver->set_box_constraints(this->get_box_constraints());
            ok = impl_->solver->solve(*this->get_operator(), b, x);
        } else {
            impl_->op.initialize(make_ref(b));

            Vector x_G;
            impl_->op.create_vector(x_G);

            if (impl_->debug) {
                std::stringstream ss;
                ss << "size: " << x.size() << "\n";
                ss << "G_size: " << x_G.size() << "\n";

                b.comm().root_print(ss.str());
            }

            if (this->has_upper_bound()) {
                auto u_G = std::make_shared<Vector>();
                impl_->op.create_vector(*u_G);
                impl_->op.select(*this->upper_bound(), *u_G);

                impl_->solver->upper_bound() = u_G;
            }

            if (this->has_lower_bound()) {
                auto l_G = std::make_shared<Vector>();
                impl_->op.create_vector(*l_G);
                impl_->op.select(*this->lower_bound(), *l_G);

                impl_->solver->lower_bound() = l_G;
            }
#if 0
            if (impl_->use_preconditioner) {
                std::cout << "BDDQPSolver: Using preconditoner!!!\n";
                auto &&interface_dofs = impl_->op.interface_dofs();
                Matrix A_GG_copy = *impl_->op.reduced_matrix();

                auto s = const_local_view_device(*impl_->boolean_selector);

                IndexSet idx;
                idx.reserve(interface_dofs.size());

                std::size_t n = interface_dofs.size();

                for (std::size_t i = 0; i < n; ++i) {
                    auto dof = interface_dofs[i];
                    if (s.get(dof) > 0) {
                        idx.push_back(i);
                    }
                }

                A_GG_copy.set_zero_rows(idx, 1);

                auto prec = std::make_shared<Factorization<Matrix, Vector>>();
                prec->update(std::make_shared<Matrix>(std::move(A_GG_copy)));

                Vector prec_residual(layout(impl_->op.righthand_side()), 0);
                prec->apply(impl_->op.righthand_side(), prec_residual);

                auto composite_op = utopia::mult<Vector>(prec, make_ref(impl_->op));
                ok = impl_->solver->solve(*composite_op, prec_residual, x_G);

            } else
#endif
            { ok = impl_->solver->solve(impl_->op, impl_->op.righthand_side(), x_G); }

            impl_->op.finalize(x_G, x);
        }

        UTOPIA_TRACE_REGION_END("BDDQPSolver::apply");
        return ok;
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::set_solver(const std::shared_ptr<MatrixFreeQPSolver> &solver) {
        impl_->solver = solver;
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::init_memory(const Layout &layout) {
        Super::init_memory(layout);

        assert(impl_->solver);

        if (layout.comm().size() == 1) {
            impl_->solver->init_memory(layout);
        } else {
            impl_->solver->init_memory(impl_->op.vector_layout());
        }
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        UTOPIA_TRACE_REGION_BEGIN("BDDQPSolver::update");

        Super::update(op);

        if (op->comm().size() == 1) {
        } else {
            if (impl_->op.selector().empty() || impl_->determine_boolean_selector) {
                determine_boolean_selector();
            }

            impl_->op.initialize(op);
        }

        init_memory(row_layout(*op));

        UTOPIA_TRACE_REGION_END("BDDQPSolver::update");
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::set_selection(const std::shared_ptr<Vector> &boolean_selector) {
        impl_->boolean_selector = boolean_selector;
        if (!boolean_selector) return;

        auto selector_view = local_view_device(*boolean_selector);

        auto &&selection = impl_->op.selector();
        if (SizeType(selection.size()) != boolean_selector->local_size()) {
            selection.resize(boolean_selector->local_size());
        }

        parallel_for(local_range_device(*boolean_selector),
                     // FIXME
                     // UTOPIA_LAMBDA(const SizeType i)
                     [&](const SizeType i) { selection[i] = selector_view.get(i) > 0; });
    }

    template <class Matrix, class Vector, int Backend>
    void BDDQPSolver<Matrix, Vector, Backend>::determine_boolean_selector() {
        UTOPIA_TRACE_REGION_BEGIN("BDDQPSolver::determine_boolean_selector");

        // std::cout << "infinity: " << impl_->infinity << "\n";

        auto selector = std::make_shared<Vector>();
        bool ok = this->get_box_constraints().determine_boolean_selector(-impl_->infinity, impl_->infinity, *selector);
        assert(ok);

        if (ok) {
            set_selection(selector);
        }

        UTOPIA_TRACE_REGION_END("BDDQPSolver::determine_boolean_selector");
    }

    template class BDDQPSolver<PetscMatrix, PetscVector, PETSC>;
}  // namespace utopia
