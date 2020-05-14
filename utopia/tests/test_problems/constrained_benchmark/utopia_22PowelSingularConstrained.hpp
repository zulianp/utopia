#ifndef UTOPIA_SOLVER_POWEL_EXTENDED_SINGULAR_CONSTRAINED
#define UTOPIA_SOLVER_POWEL_EXTENDED_SINGULAR_CONSTRAINED

#include <cassert>
#include <vector>
#include "utopia_Function.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class ExtendedPowell22Constrained final : public ConstrainedTestFunction<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Matrix>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        ExtendedPowell22Constrained() {
            auto v_layout = serial_layout(dim());

            Vector ub, lb;
            ub.zeros(v_layout);
            lb.zeros(v_layout);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);

                lb.set(0, 0.1);
                lb.set(1, -20.0);
                lb.set(2, -1.0);
                lb.set(3, -1.0);

                ub.set(0, 100.0);
                ub.set(1, 20.0);
                ub.set(2, 1.0);
                ub.set(3, 50.0);
            }

            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));
        }

        std::string name() const override { return "Extended Powell singular, bound constrained"; }

        SizeType dim() const override { return unconstrained_.dim(); }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            return unconstrained_.value(x, result);
        }

        bool gradient(const Vector &x, Vector &g) const override { return unconstrained_.gradient(x, g); }

        bool hessian(const Vector &x, Matrix &H) const override { return unconstrained_.hessian(x, H); }

        Vector initial_guess() const override { return unconstrained_.initial_guess(); }

        const Vector &exact_sol() const override { return unconstrained_.exact_sol(); }

        Scalar min_function_value() const override { return 0.1878196e-3; }

    private:
        ExtendedPowell22<Matrix, Vector> unconstrained_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_POWEL_EXTENDED_SINGULAR_CONSTRAINED
