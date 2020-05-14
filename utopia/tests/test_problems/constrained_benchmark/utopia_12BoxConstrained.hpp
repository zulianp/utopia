#ifndef UTOPIA_BOX_12_CONSTRAINED
#define UTOPIA_BOX_12_CONSTRAINED

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Box12Constrained final : public ConstrainedTestFunction<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Matrix>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        Box12Constrained() {
            auto v_layout = serial_layout(dim());

            Vector ub, lb;
            ub.zeros(v_layout);
            lb.zeros(v_layout);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);

                lb.set(0, 0.0);
                lb.set(1, 5.0);
                lb.set(2, 0.0);

                ub.set(0, 2.0);
                ub.set(1, 9.5);
                ub.set(2, 20.0);
            }

            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));
        }

        std::string name() const override { return "Box three-dimensional, bound constrained"; }

        SizeType dim() const override { return unconstrained_.dim(); }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            return unconstrained_.value(x, result);
        }

        bool gradient(const Vector &x, Vector &g) const override { return unconstrained_.gradient(x, g); }

        bool hessian(const Vector &x, Matrix &H) const override { return unconstrained_.hessian(x, H); }

        Vector initial_guess() const override { return unconstrained_.initial_guess(); }

        const Vector &exact_sol() const override { return unconstrained_.exact_sol(); }

        Scalar min_function_value() const override { return 0.30998153e-5; }

    private:
        Box12<Matrix, Vector> unconstrained_;
    };
}  // namespace utopia

#endif  // UTOPIA_BOX_12_CONSTRAINED