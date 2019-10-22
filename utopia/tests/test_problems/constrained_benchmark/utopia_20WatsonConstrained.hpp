#ifndef UTOPIA_SOLVER_WATSON_20_CONSTRAINED_HPP
#define UTOPIA_SOLVER_WATSON_20_CONSTRAINED_HPP

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Watson20Constrained final: public ConstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Watson20Constrained()
        {
            assert(mpi_world_size() == 1 && "does not work for parallel matrices");

            Vector ub, lb;
            ub = zeros(9);
            lb = zeros(9);

            {
                const Write<Vector> write2(ub);
                const Write<Vector> write1(lb);

                lb.set(0, -0.00001);
                lb.set(1, 0.0);
                lb.set(2, 0.0);
                lb.set(3, 0.0);
                lb.set(4, 0.0);
                lb.set(5, -3.0);
                lb.set(6, 0.0);
                lb.set(7, -3.0);
                lb.set(8, 0.0);

                ub.set(0, 0.00001);
                ub.set(1, 0.9);
                ub.set(2, 0.1);
                ub.set(3, 1.0);
                ub.set(4, 1.0);
                ub.set(5, 0.0);
                ub.set(6, 4.0);
                ub.set(7, 0.0);
                ub.set(8, 2.0);                
            }

            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));
        }

        std::string name() const override
        {
            return "Watson, bound constrained_1";
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
            return 0.37401397e-1; 
        }

    private:
        Watson20<Matrix, Vector> unconstrained_; 
    };

}

#endif //UTOPIA_SOLVER_WATSON_20_CONSTRAINED_HPP
