#ifndef UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP
#define UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP

#include <vector>
#include <assert.h>
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
    class TestFunction2D_1 : public Function<Matrix, Vector> 
    {
    public:

        TestFunction2D_1() { };

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

        bool hessian(const Vector &point, Matrix &result) const override {

            result = zeros(2, 2);

            const Write<Matrix> write(result);

            result.set(0, 0, 4.0);
            result.set(1, 1, 8.0);
            return true;
        }
    };


    /**
     * @brief      Rosenbrock 2D banana function. \n 
     *             The floor of the valley follows approximately the parabola \f$ y = x^2 + 1/200 \f$.   
     *             The covariance matrix is not positive-definite. On the dashed line it is singular. 
     *             Stepping method tend to perform at least as well as gradient methods for this function.
     *
     */
    template<class Matrix, class Vector>
    class Rosenbrock : public Function<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        Rosenbrock() {
            //FIXME find a way to implement generic Rosenbrock in parallel
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result = 1 + 100.0 * pow(x * x - y , 2.0) + pow(x - 1 , 2.0);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override 
        {
            assert(point.size().get(0) == 2);
            result = zeros(2);

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result.set(0, (400.0 * x * x * x - 400 * x * y + 2.0 * x - 2.0));
            result.set(1, 200.0 * (y - x * x));
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override 
        {
            assert(point.size().get(0) == 2);


            result = zeros(2, 2);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar mixed = -400.0 * x;

            result.set(0, 0, 1200 * x * x - 400 * y + 2);
            result.set(0, 1, mixed);
            result.set(1, 0, mixed);
            result.set(1, 1, 200.0);
            return true;
        }
    };
}

#endif //UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP
