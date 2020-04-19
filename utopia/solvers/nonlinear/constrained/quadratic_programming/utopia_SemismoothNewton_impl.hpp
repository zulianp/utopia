#ifndef UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP
#define UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP

#include "utopia_SemismoothNewton.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class SemismoothNewton<Matrix, Vector, Backend>::Buffers {
    public:
        void init(const Layout &layout) {
            correction.zeros(layout);
            residual.zeros(layout);

            d_lb.zeros(layout);
            d_ub.zeros(layout);

            lambda_lb.zeros(layout);
            lambda_ub.zeros(layout);

            active.zeros(layout);
        }

        Vector correction, residual;
        Vector d_lb, d_ub;
        Vector lambda_lb, lambda_ub;
        Vector active;

        Matrix mat;
    };

    template <class Matrix, class Vector, int Backend>
    SemismoothNewton<Matrix, Vector, Backend>::SemismoothNewton(const std::shared_ptr<LinearSolver> &linear_solver)
        : linear_solver_(linear_solver), buffers_(utopia::make_unique<Buffers>()) {}

    template <class Matrix, class Vector, int Backend>
    SemismoothNewton<Matrix, Vector, Backend>::~SemismoothNewton() {}

    template <class Matrix, class Vector, int Backend>
    SemismoothNewton<Matrix, Vector, Backend> *SemismoothNewton<Matrix, Vector, Backend>::clone() const {
        auto ptr = utopia::make_unique<SemismoothNewton>(std::shared_ptr<LinearSolver>(linear_solver_->clone()));
        if (this->has_bound()) {
            ptr->set_box_constraints(this->get_box_constraints());
        }
        return ptr.release();
    }

    template <class Matrix, class Vector, int Backend>
    void SemismoothNewton<Matrix, Vector, Backend>::read(Input &in) {
        QPSolver<Matrix, Vector>::read(in);
        if (linear_solver_) {
            in.get("linear_solver", *linear_solver_);
        }
    }

    template <class Matrix, class Vector, int Backend>
    void SemismoothNewton<Matrix, Vector, Backend>::print_usage(std::ostream &os) const {
        Super::print_usage(os);
    }

    template <class Matrix, class Vector, int Backend>
    bool SemismoothNewton<Matrix, Vector, Backend>::smooth(const Vector &, Vector &) {
        assert(false && "IMPLEMENT ME");
        return false;
    }

    template <class Matrix, class Vector, int Backend>
    bool SemismoothNewton<Matrix, Vector, Backend>::apply(const Vector &b, Vector &x) {
        const auto &A = *this->get_operator();

        assert(A.comm().size() == b.comm().size());
        assert(A.rows() == b.size());
        assert(x.empty() || (x.size() == b.size() && b.comm().size() == x.comm().size()));
        assert(this->has_bound());

        this->fill_empty_bounds(layout(b));

        auto &&lb = this->get_lower_bound();
        auto &&ub = this->get_upper_bound();

        auto &correction = buffers_->correction;
        auto &residual = buffers_->residual;

        auto &d_lb = buffers_->d_lb;
        auto &d_ub = buffers_->d_ub;

        auto &lambda_lb = buffers_->lambda_lb;
        auto &lambda_ub = buffers_->lambda_ub;

        auto &mat = buffers_->mat;

        auto &active = buffers_->active;

        lambda_lb.set(0.0);
        lambda_ub.set(0.0);

        bool converged = false;

        auto r = local_range_device(b);

        if (this->verbose()) {
            this->init_solver("SemismoothNewton comm.size = " + std::to_string(b.comm().size()),
                              {" it. ", "      || x_k - x_{k-1} ||"});
        }

        SizeType it = 0;

        while (!converged) {
            d_lb = lambda_lb + x - lb;
            d_ub = lambda_ub + x - ub;

            auto d_lb_view = const_local_view_device(d_lb);
            auto d_ub_view = const_local_view_device(d_ub);
            auto active_view = local_view_device(active);

            SizeType changed = 0;
            parallel_reduce(r,
                            UTOPIA_LAMBDA(const SizeType &i)->SizeType {
                                const bool prev_active = active_view.get(i);

                                const Scalar d_lb_i = d_lb_view.get(i);
                                const Scalar d_ub_i = d_ub_view.get(i);

                                const bool active_lb = d_lb_i <= 0.0;  // FIXME use tol
                                const bool active_ub = d_ub_i >= 0.0;  // FIXME use tol
                                const bool active_i = active_lb || active_ub;

                                active_view.set(i, active_i);
                                return prev_active != active_i;
                            },
                            changed);

            residual = A * x;
            residual = b - residual;

            residual -= e_mul(active, residual);

            // maybe there is a more efficent way than copy everything but this is probably ok
            mat.same_nnz_pattern_copy(A);
            set_zero_rows(mat, active, 1.0);

            linear_solver_->solve(mat, residual, correction);

            x += correction;

            Scalar norm_c = norm2(correction);
            ++it;

            if (this->verbose()) {
                PrintInfo::print_iter_status(it, {norm_c});
            }

            if (it > 1 && !changed) {
                // we converged
                converged = true;
            } else {
                converged = this->check_convergence(it, 1, 1, norm_c);
            }
        }

        return converged;
    }

    template <class Matrix, class Vector, int Backend>
    void SemismoothNewton<Matrix, Vector, Backend>::init_memory(const Layout &layout) {
        Super::init_memory(layout);
        if (linear_solver_) {
            linear_solver_->init_memory(layout);
        }

        buffers_->init(layout);
    }

    template <class Matrix, class Vector, int Backend>
    void SemismoothNewton<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        Super::update(op);
        init_memory(row_layout(*op));

        // copy the matrix
        buffers_->mat = *op;
    }

}  // namespace utopia

#endif  // UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP