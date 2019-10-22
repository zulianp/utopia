#ifndef UTOPIA_SOLVER_WOODS_CONSTRAINED_HPP
#define UTOPIA_SOLVER_WOODS_CONSTRAINED_HPP

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Woods14Constrained final: public ConstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Woods14Constrained()
        {
            assert(mpi_world_size() == 1 && "does not work for parallel matrices");

            Vector ub, lb;
            ub = zeros(4);
            lb = zeros(4);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);
                
                lb.set(0, -100.0);
                lb.set(1, -100.0);
                lb.set(2, -100.0);
                lb.set(3, -100.0);

                ub.set(0, 0.0);
                ub.set(1, 10);
                ub.set(2, 100);
                ub.set(3, 100);
            }

            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));
        }

        std::string name() const override
        {
            return "Woods, bound constrained.";
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
            return 0.15567008e1;
        }

    private:
        Woods14<Matrix, Vector> unconstrained_; 

    };



}

#endif //UTOPIA_SOLVER_WOODS_CONSTRAINED_HPP
