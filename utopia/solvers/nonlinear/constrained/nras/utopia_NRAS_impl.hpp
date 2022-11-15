#ifndef UTOPIA_NRAS_IMPL_HPP
#define UTOPIA_NRAS_IMPL_HPP

// https://petsc.org/release/src/dm/tutorials/ex19.c.html

#include "utopia_NRAS.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class NRAS<Matrix, Vector>::Impl {
    public:
        std::shared_ptr<Function> global_function;
        std::shared_ptr<Function> local_function;

        std::shared_ptr<Operator<Vector>> cutoff;
        std::shared_ptr<Operator<Vector>> projection;
    };

    template <class Matrix, class Vector>
    NRAS<Matrix, Vector>::NRAS() : impl_(utopia::make_unique<Impl>()) {}

    template <class Matrix, class Vector>
    NRAS<Matrix, Vector>::~NRAS() {}

    template <class Matrix, class Vector>
    void NRAS<Matrix, Vector>::set_local_function(const std::shared_ptr<Function> &lfun) {
        impl_->local_function = lfun;
    }

    template <class Matrix, class Vector>
    void NRAS<Matrix, Vector>::set_global_function(const std::shared_ptr<Function> &gfun) {
        impl_->global_function = gfun;
    }

    template <class Matrix, class Vector>
    void NRAS<Matrix, Vector>::set_cutoff_operator(const std::shared_ptr<Operator<Vector>> &op) {
        impl_->cutoff = op;
    }

    template <class Matrix, class Vector>
    void NRAS<Matrix, Vector>::set_projection_operator(const std::shared_ptr<Operator<Vector>> &op) {
        impl_->projection = op;
    }

    template <class Matrix, class Vector>
    bool NRAS<Matrix, Vector>::solve(Vector &x) {
        Vector local_grad_x;
        Vector local_view_x;
        Vector replicated_x;

        auto &&cutoff = *impl_->cutoff;
        auto &&projection = *impl_->projection;
        auto &&lfun = *impl_->local_function;
        auto &&gfun = *impl_->global_function;

        cutoff.apply(x, replicated_x);

        replicated_x.create_local_vector(local_view_x);

        disp(local_view_x);

        lfun.gradient(local_view_x, local_grad_x);

        for (int r = 0; r < x.comm().size(); ++r) {
            if (r == x.comm().rank()) {
                disp(local_grad_x);
            }

            x.comm().barrier();
        }

        replicated_x.restore_local_vector(local_view_x);
        return false;
    }
}  // namespace utopia

#endif
