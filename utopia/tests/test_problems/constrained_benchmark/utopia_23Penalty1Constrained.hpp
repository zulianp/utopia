#ifndef UTOPIA_SOLVER_PENALTY1_23_CONSTRAINED_HPP
#define UTOPIA_SOLVER_PENALTY1_23_CONSTRAINED_HPP

#include <vector>
#include <assert.h>
#include "utopia_UnconstrainedBenchmark.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class PenaltyI23Constrained final: public ConstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        PenaltyI23Constrained()
        {
            auto v_layout = serial_layout(dim());

            Vector ub, lb;
            ub.zeros(v_layout);
            lb.zeros(v_layout);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);

                lb.set(0, 0.0);
                lb.set(1, 1.0);
                lb.set(2, 0.0);
                lb.set(3, 0.0);
                lb.set(4, 0.0);
                lb.set(5, 1.0);
                lb.set(6, 0.0);
                lb.set(7, 0.0);
                lb.set(8, 0.0);
                lb.set(9, 1.0);

                ub.set(0, 100.0);
                ub.set(1, 100.0);
                ub.set(2, 100.0);
                ub.set(3, 100.0);
                ub.set(4, 100.0);
                ub.set(5, 100.0);
                ub.set(6, 100.0);
                ub.set(7, 100.0);
                ub.set(8, 100.0);
                ub.set(9, 100.0);
            }


            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));

        }

        std::string name() const override
        {
            return "Penalty I, bound constrained";
        }

        SizeType dim() const override
        {
            return unconstrained_.dim();
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            return unconstrained_.value(x, result);
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            return unconstrained_.gradient(x, g);
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            return unconstrained_.hessian(x, H);
        }

        Vector initial_guess() const override
        {
            return unconstrained_.initial_guess();
        }

        const Vector & exact_sol() const override
        {
            return unconstrained_.exact_sol();
        }

        Scalar min_function_value() const override
        {
            return 0.75625699e1;
        }

    private:
        PenaltyI23<Matrix, Vector> unconstrained_;

    };

}

#endif //UTOPIA_SOLVER_PENALTY1_23_CONSTRAINED_HPP
