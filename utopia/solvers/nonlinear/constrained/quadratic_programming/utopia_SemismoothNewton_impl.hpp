#ifndef UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP
#define UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP

#include "utopia_SemismoothNewton.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class SemismoothNewton<Matrix, Vector, Backend>::Buffers {
    public:
        void init(const Layout &layout) {
            x_old.zeros(layout);
            d_lb.zeros(layout);
            d_ub.zeros(layout);
            lambda_lb.zeros(layout);
            lambda_ub.zeros(layout);
        }

        Vector x_old;
        Vector d_lb, d_ub;
        Vector lambda_lb, lambda_ub;
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

        auto &x_old = buffers_->x_old;
        auto &d_lb = buffers_->d_lb;
        auto &d_ub = buffers_->d_ub;
        auto &lambda_lb = buffers_->lambda_lb;
        auto &lambda_ub = buffers_->lambda_ub;

        lambda_lb.set(0.0);
        lambda_ub.set(0.0);

        bool converged = false;

        auto r = local_range_device(b);

        while (!converged) {
            d_lb = lambda_lb + x - lb;
            d_ub = lambda_ub + x - ub;

            auto d_lb_view = const_local_view_device(d_lb);
            auto d_ub_view = const_local_view_device(d_ub);

            parallel_for(r,
                         UTOPIA_LAMBDA(const SizeType &i){

                         });
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
    }

}  // namespace utopia

#endif  // UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP