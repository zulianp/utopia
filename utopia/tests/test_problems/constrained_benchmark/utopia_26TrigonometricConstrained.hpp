#ifndef UTOPIA_SOLVER_TRIGONOMETRIC_26_CONSTRAINED_HPP
#define UTOPIA_SOLVER_TRIGONOMETRIC_26_CONSTRAINED_HPP

#include <cassert>
#include <vector>
#include "utopia_UnconstrainedBenchmark.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Trigonometric26Constrained final : public ConstrainedTestFunction<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Matrix>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        Trigonometric26Constrained() {
            auto v_layout = serial_layout(dim());

            Vector ub, lb;
            ub.zeros(v_layout);
            lb.zeros(v_layout);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);

                lb.set(0, 0.0);
                lb.set(1, 10.0);
                lb.set(2, 20.0);
                lb.set(3, 30.0);
                lb.set(4, 40.0);
                lb.set(5, 50.0);
                lb.set(6, 60.0);
                lb.set(7, 70.0);
                lb.set(8, 80.0);
                lb.set(9, 90.0);

                ub.set(0, 10.0);
                ub.set(1, 20.0);
                ub.set(2, 30.0);
                ub.set(3, 40.0);
                ub.set(4, 50.0);
                ub.set(5, 60.0);
                ub.set(6, 70.0);
                ub.set(7, 80.0);
                ub.set(8, 90.0);
                ub.set(9, 100.0);
            }

            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));
        }

        std::string name() const override { return "Trigonometric, bound constrained"; }

        SizeType dim() const override { return unconstrained_.dim(); }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            return unconstrained_.value(x, result);
        }

        bool gradient(const Vector &x, Vector &g) const override { return unconstrained_.gradient(x, g); }

        bool hessian(const Vector &x, Matrix &H) const override { return unconstrained_.hessian(x, H); }

        Vector initial_guess() const override { return unconstrained_.initial_guess(); }

        const Vector &exact_sol() const override { return unconstrained_.exact_sol(); }

    private:
        Trigonometric26<Matrix, Vector> unconstrained_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_TRIGONOMETRIC_26_CONSTRAINED_HPP
