#ifndef UTOPIA_SHIFTED_PENALTY_QP_SOLVER_IMPL_HPP
#define UTOPIA_SHIFTED_PENALTY_QP_SOLVER_IMPL_HPP

#include "utopia_ShiftedPenaltyQPSolver.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

#include "utopia_ShiftedPenalty.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::Impl {
    public:
        void init(const Layout &layout) {
            if (linear_solver) {
                linear_solver->init_memory(layout);
            }
        }

        Impl() : penalty(utopia::make_unique<ShiftedPenalty<Matrix, Vector>>()) {}

        bool debug{false};

        std::shared_ptr<LinearSolver> linear_solver;
        std::shared_ptr<Vector> boolean_selector;
        std::shared_ptr<Matrix> scaling_matrix;
        std::unique_ptr<Penalty<Matrix, Vector>> penalty;
        std::shared_ptr<Transformation> transform;
    };

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::set_scaling_matrix(
        const std::shared_ptr<Matrix> &scaling_matrix) {
        this->impl_->scaling_matrix = scaling_matrix;
    }

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
        impl_->penalty->read(in);

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

        // this->verbose(true);

        if (this->verbose()) {
            this->init_solver("ShiftedPenaltyQPSolver comm.size = " + std::to_string(b.comm().size()),
                              {" it. ",
                               "      || x_k - x_{k-1} ||, "
                               "      || g_{k-1} ||, "});
        }

        Vector c, g, d;
        Matrix H, H_penalty;

        H = *this->get_operator();
        d = diag(H);

        c.zeros(layout(b));
        g.zeros(layout(b));

        auto H_unconstrained_ptr = this->get_operator();
        auto linear_solver = impl_->linear_solver;

        assert(impl_->penalty);
        auto &&penalty = impl_->penalty;
        penalty->reset();
        penalty->set_box_constraints(make_ref(this->get_box_constraints()));
        penalty->set_selection(impl_->boolean_selector);
        penalty->set_scaling_matrix(impl_->scaling_matrix);
        penalty->set_transform(impl_->transform);

        bool converged = false;

        const int max_it = this->max_it();
        for (int it = 1; it <= max_it; it++) {
            g = *H_unconstrained_ptr * x - b;

            // Set buffers to zero
            if (!c.empty()) c.set(0);
            if (!H_penalty.empty()) H_penalty *= 0;

            penalty->update(x);
            penalty->gradient(x, c);
            penalty->hessian(x, H_penalty);

            g += c;
            g = -g;

            // TODO Find a way to be more efficient
            H = *H_unconstrained_ptr;
            H += H_penalty;

            // disp(H_penalty);

            c.set(0);
            linear_solver->solve(H, g, c);
            x += c;

            const Scalar norm_c = norm2(c);
            if (this->verbose()) {
                const Scalar norm_g = norm2(g);
                PrintInfo::print_iter_status(it, {norm_c, norm_g});
            }

            converged = this->check_convergence(it, 1, 1, norm_c);
            if (converged) break;
        }

        return converged;
    }

    template <class Matrix, class Vector, int Backend>
    void ShiftedPenaltyQPSolver<Matrix, Vector, Backend>::set_transform(const std::shared_ptr<Transformation> &t) {
        impl_->transform = t;
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
        impl_->boolean_selector = boolean_selector;
    }

}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_QP_SOLVER_IMPL_HPP