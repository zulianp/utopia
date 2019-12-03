#ifndef UTOPIA_SOLVER_TESTFUNCTION2D_HPP
#define UTOPIA_SOLVER_TESTFUNCTION2D_HPP

#include <vector>
#include <assert.h>
#include <cmath>
#include "utopia_Function.hpp"



namespace utopia
{

    /**
     * @brief      Example of the nonlinear function. Used to test nonlinear solvers.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class QPTestFunction_2D final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:

        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        QPTestFunction_2D() 
        {
            x_init_ = zeros(2);
            x_exact_ = zeros(2);
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            const Read<Vector> read(point);

            result = 4 * ((3.0 - 0.5 * point.get(0)) * (3.0 - 0.5 * point.get(0)) +
                        (point.get(1) + 7.0) * (point.get(1) + 7.0));

            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {

            result = zeros(2);

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            result.set(0, 4.0 * (0.5 * point.get(0) - 3.0));
            result.set(1, 8.0 * (point.get(1) + 7.0));
            return true;
        }

        bool hessian(const Vector &/*point*/, Matrix &result) const override {

            result = zeros(2, 2);

            const Write<Matrix> write(result);

            result.set(0, 0, 4.0);
            result.set(1, 1, 8.0);
            return true;
        }


        Vector initial_guess() const override
        {
            return x_init_;
        }

        const Vector & exact_sol() const override
        {
            return x_exact_;
        }

        Scalar min_function_value() const override
        {
            return 1; // TBD
        }

        std::string name() const override
        {
            return "QPTestFunction_2D";
        }

        SizeType dim() const override
        {
            return 2;
        }

        bool exact_sol_known() const
        {
            return false;
        }        


    private:
        Vector x_init_;
        Vector x_exact_;        
    };

}

#endif //UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP
