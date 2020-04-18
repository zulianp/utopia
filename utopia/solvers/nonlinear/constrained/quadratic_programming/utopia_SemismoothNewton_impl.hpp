#ifndef UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP
#define UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP
#include "utopia_SemismoothNewton.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    SemismoothNewton<Matrix, Vector, Backend>::SemismoothNewton(const std::shared_ptr<LinearSolver> &linear_solver)
        : linear_solver_(linear_solver) {}

    template <class Matrix, class Vector, int Backend>
    SemismoothNewton<Matrix, Vector, Backend> *SemismoothNewton<Matrix, Vector, Backend>::clone() const {
        auto ptr = utopia::make_unique<SemismoothNewton>(std::shared_ptr<LinearSolver>(linear_solver_->clone()));
        // if (this->has_bound()) {
        //     ptr->set_box_constraints(this->get_box_constraints());
        // }
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

        return true;
    }

    template <class Matrix, class Vector, int Backend>
    void SemismoothNewton<Matrix, Vector, Backend>::init_memory(const Layout &layout) {
        Super::init_memory(layout);
        if (linear_solver_) {
            linear_solver_->init_memory(layout);
        }
    }

    template <class Matrix, class Vector, int Backend>
    void SemismoothNewton<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        Super::update(op);

        if (linear_solver_) {
            linear_solver_->init_memory(row_layout(*op));
        }
    }

}  // namespace utopia

#endif  // UTOPIA_SEMISMOOTH_NEWTON_IMPL_HPP