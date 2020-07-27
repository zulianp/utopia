#ifndef UTOPIA_TRILINOS_DIFF_CONTROLLER_HPP
#define UTOPIA_TRILINOS_DIFF_CONTROLLER_HPP

#include "utopia_trilinos_Base.hpp"

// FIXME The default implementation does not work with cuda
// Provide cuda compatible implementation here
#ifdef KOKKOS_ENABLE_CUDA

#include "utopia_DiffController.hpp"
#include "utopia_trilinos_ForwardDeclarations.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class DiffController<Matrix, Vector, TRILINOS> : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        DiffController(const Scalar spacing = 1e-5) : spacing_(spacing) {}

        inline void hessian_from_grad(const bool val) { hessian_from_grad_ = val; }

        void read(Input &in) override {
            in.get("spacing", spacing_);
            in.get("hessian_from_grad", hessian_from_grad_);
        }

        template <class Fun>
        bool check(const Fun &fun, const Vector &x, const Vector &g, const Matrix &H) const {
            return check_grad(fun, x, g) && check_hessian(fun, x, H);
        }

        template <class Fun>
        bool check_grad(const Fun &fun, const Vector &x, const Vector &g) const {
            // IMPLEMENT ME
            return true;
        }

        template <class Fun>
        bool check_hessian(const Fun &fun, const Vector &x, const Matrix &H) const {
            // IMPLEMENT ME
            return true;
        }

        void spacing(const Scalar &s) { spacing_ = s; }

    private:
        Scalar spacing_;
        bool hessian_from_grad_{true};
    };
}  // namespace utopia

#endif  // KOKKOS_ENABLE_CUDA
#endif  // UTOPIA_TRILINOS_DIFF_CONTROLLER_HPP
