#ifndef UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
#define UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP

#include <memory>
#include "utopia_QPSolver.hpp"
#include "utopia_SemismoothNewton.hpp"
#include "utopia_petsc.hpp"
#include "utopia_petsc_impl.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> final : public QPSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::QPSolver<Matrix, Vector> Super;
        using BoxConstraints = utopia::BoxConstraints<Vector>;

    public:
        SemismoothNewton(
            const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<Factorization<Matrix, Vector>>());
        ~SemismoothNewton() override;
        SemismoothNewton *clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void read(Input &in) override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
