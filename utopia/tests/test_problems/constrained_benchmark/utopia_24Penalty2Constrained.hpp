#ifndef UTOPIA_SOLVER_PENALTY_II_CONSTRAINED
#define UTOPIA_SOLVER_PENALTY_II_CONSTRAINED

#include <cassert>
#include <vector>
#include "utopia_UnconstrainedBenchmark.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class PenaltyII24Constrained final : public ConstrainedTestFunction<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Matrix>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        PenaltyII24Constrained() {
            auto v_layout = serial_layout(dim());

            Vector ub, lb;
            ub.zeros(v_layout);
            lb.zeros(v_layout);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);

                lb.set(0, -10.0);
                lb.set(1, 0.1);
                lb.set(2, 0.0);
                lb.set(3, 0.05);
                lb.set(4, 0.0);
                lb.set(5, -10.0);
                lb.set(6, 0.0);
                lb.set(7, 0.2);
                lb.set(8, 0.0);
                lb.set(9, 0.0);

                ub.set(0, 50.0);
                ub.set(1, 50.0);
                ub.set(2, 50.0);
                ub.set(3, 50.0);
                ub.set(4, 50.0);
                ub.set(5, 50.0);
                ub.set(6, 50.0);
                ub.set(7, 50.0);
                ub.set(8, 50.0);
                ub.set(9, 0.5);
            }

            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));
        }

        std::string name() const override { return "Penalty II, bound constrained"; }

        SizeType dim() const override { return unconstrained_.dim(); }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            return unconstrained_.value(x, result);
        }

        bool gradient(const Vector &x, Vector &g) const override { return unconstrained_.gradient(x, g); }

        bool hessian(const Vector &x, Matrix &H) const override { return unconstrained_.hessian(x, H); }

        Vector initial_guess() const override { return unconstrained_.initial_guess(); }

        const Vector &exact_sol() const override { return unconstrained_.exact_sol(); }

        Scalar min_function_value() const override { return 0.29442600e-3; }

    private:
        PenaltyII24<Matrix, Vector> unconstrained_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_PENALTY_II_CONSTRAINED
