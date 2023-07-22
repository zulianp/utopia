#ifndef UTOPIA_SHIFTED_PENALTY_QP_SOLVER_IMPL_HPP
#define UTOPIA_SHIFTED_PENALTY_QP_SOLVER_IMPL_HPP

#include "utopia_ShiftedPenaltyQPSolver.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::Impl {
    public:
        void init(const Layout &layout) {
            if (linear_solver) {
                linear_solver->init_memory(layout);
            }
        }

        std::shared_ptr<LinearSolver> linear_solver;
        bool debug{false};
        Scalar penalty_param{1e6};
    };

    template <class Matrix, class Vector, int Backend>
    ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::ShiftedPenaltyQPSolver(std::shared_ptr<LinearSolver> linear_solver)
        : impl_(utopia::make_unique<Impl>()) {
        set_linear_solver(linear_solver);
    }

    template <class Matrix, class Vector, int Backend>
    ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::ShiftedPenaltyQPSolver() : impl_(utopia::make_unique<Impl>()) {
        set_linear_solver(std::make_shared<OmniLinearSolver<Matrix, Vector>>());
    }

    template <class Matrix, class Vector, int Backend>
    ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::~ShiftedPenaltyQPSolver() = default;

    template <class Matrix, class Vector, int Backend>
    ShiftedPenaltyQPSolver<Matrix, Vector, Backend> *ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::clone() const {
        auto ptr =
            utopia::make_unique<ShiftedPenaltyQPSolver>(std::shared_ptr<LinearSolver>(impl_->linear_solver->clone()));
        if (this->has_bound()) {
            ptr->set_box_constraints(this->get_box_constraints());
        }
        return ptr.release();
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::read(Input &in) {
        QPSolver<Matrix, Vector>::read(in);
        in.get("debug", impl_->debug);
        in.get("penalty_param", impl_->penalty_param);

        if (impl_->linear_solver) {
            in.get("linear_solver", *impl_->linear_solver);
        }
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::print_usage(std::ostream &os) const {
        Super::print_usage(os);
    }

    template <class Matrix, class Vector, int Backend>
    bool ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::smooth(const Vector &, Vector &) {
        assert(false && "IMPLEMENT ME");
        return false;
    }

    template <class Matrix, class Vector, int Backend>
    bool ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::apply(const Vector &b, Vector &x) {
        assert(this->get_operator()->comm().size() == b.comm().size());
        assert(this->get_operator()->rows() == b.size());
        assert(x.empty() || (x.size() == b.size() && b.comm().size() == x.comm().size()));
        assert(this->has_upper_bound());

        this->verbose(true);

        if (this->verbose()) {
            this->init_solver("ShiftedPenaltyQPSolver comm.size = " + std::to_string(b.comm().size()),
                              {" it. ", "      || x_k - x_{k-1} ||"});
        }

        Vector d, g, diag_B, c;
        Matrix A;

        A = *this->get_operator();
        diag_B.zeros(layout(b));
        c.zeros(layout(b));

        auto ub_ptr = this->upper_bound();
        auto a_ptr = this->get_operator();
        auto linear_solver = impl_->linear_solver;

        const Scalar penalty_param = impl_->penalty_param;

        bool converged = false;

        const int max_it = this->max_it();
        for (int it = 1; it <= max_it; it++) {
            g = *a_ptr * x - b;
            d = *ub_ptr - x;

            A.same_nnz_pattern_copy(*this->get_operator());

            {
                auto d_view = const_local_view_device(d);
                auto g_view = local_view_device(g);
                auto diag_B_view = local_view_device(diag_B);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar gi = g_view.get(i);
                        const Scalar active = di <= 0;
                        const Scalar g_active = active * (gi - penalty_param * di);
                        const Scalar gg = -(gi + g_active);
                        g_view.set(i, gg);
                        diag_B_view.set(i, active * penalty_param);
                    });
            }

            A.shift_diag(diag_B);
            c.set(0);
            linear_solver->solve(A, g, c);
            x += c;

            Scalar norm_c = norm2(c);

            if (this->verbose()) {
                PrintInfo::print_iter_status(it, {norm_c});
            }

            converged = this->check_convergence(it, 1, 1, norm_c);
            if (converged) break;
        }

        return converged;
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::set_linear_solver(
        const std::shared_ptr<LinearSolver> &linear_solver) {
        impl_->linear_solver = linear_solver;
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::init_memory(const Layout &layout) {
        Super::init_memory(layout);
        impl_->init(layout);
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        Super::update(op);
        init_memory(row_layout(*op));
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::set_selection(
        const std::shared_ptr<Vector> &boolean_selector) {
        // impl_->boolean_selector = boolean_selector;
    }

}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_QP_SOLVER_IMPL_HPP