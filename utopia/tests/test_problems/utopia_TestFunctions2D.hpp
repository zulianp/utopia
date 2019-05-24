#ifndef UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP
#define UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP

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


    template<class Matrix, class Vector>
    class SmallSingularExample : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        SmallSingularExample()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            assert(x.size().get(0) == 2);
            Vector g = values(2, 0.0);
            gradient(x, g);
            result = 0.5 * norm2(g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(x.size().get(0) == 2);
            g = zeros(2);

            auto pi = std::acos(-1.0);

            const Read<Vector> read(x);
            const Write<Vector> write(g);

            const Scalar x1 = x.get(0);
            const Scalar x2 = x.get(1);

            g.set(0, x1*x1 - x2 + 1.0);
            g.set(1, x1 - std::cos(pi/2*x2));
            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            assert(x.size().get(0) == 2);

            H = zeros(2, 2);

            const Read<Vector> read(x);
            const Write<Matrix> write(H);

            const Scalar x1 = x.get(0);
            const Scalar x2 = x.get(1);

            auto pi = std::acos(-1.0);

            H.set(0, 0, 2.0 * x1);
            H.set(0, 1, -1.0);
            H.set(1, 0, 1.0);
            H.set(1, 1, pi/2. * std::sin(pi/2.*x2));
            return true;
        }
    };


    template<class Matrix, class Vector>
    class ContinuousStirredReactor : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        ContinuousStirredReactor()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            assert(x.size().get(0) == 4);
            Vector g = values(2, 0.0);
            gradient(x, g);
            result = 0.5 * norm2(g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(x.size().get(0) == 4);
            g = zeros(4);

            {

                const Read<Vector> read(x);
                const Write<Vector> write(g);


                const Scalar x1 = x.get(0);
                const Scalar x2 = x.get(1);
                const Scalar x3 = x.get(2);
                const Scalar x4 = x.get(3);

                g.set(0, 1.0-x1-(100.0*x1*x2));
                g.set(1, 2.0-x2-(100.0*x1*x2)-(100.0*x2*x3)); 
                g.set(2, -x3+(100.0*x1*x2)-(100*x2*x3));
                g.set(3, -x4+(100.0*x2*x3));
            }

            g = -1.0*g; 

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            assert(x.size().get(0) == 4);

            H = zeros(4,4);

            {
                const Read<Vector> read(x);
                const Write<Matrix> write(H);

                const Scalar x1 = x.get(0);
                const Scalar x2 = x.get(1);
                const Scalar x3 = x.get(2);
                const Scalar x4 = x.get(3);

                H.set(0, 0, -1.0-(100.0*x2));
                H.set(0, 1, -100.0*x1);
                H.set(0, 2, 0.0);
                H.set(0, 3, 0.0);

                H.set(1, 0, -100.0*x2);
                H.set(1, 1, -1.0-(100.0*x1)-(100.*x3));
                H.set(1, 2, -100.0*x2);
                H.set(1, 3, 0.0);

                H.set(2, 0, 100.0*x2);
                H.set(2, 1, 100.0*(x1-x3));
                H.set(2, 2, -1.0-(100.0*x2));
                H.set(2, 3, 0.0);            

                H.set(3, 0, 0.0);
                H.set(3, 1, 100.0*x3);
                H.set(3, 2, 100.0*x2);
                H.set(3, 3, -1.0);
            }

            H = -1.0*H; 

            return true;
        }


        void get_initial_guess(Vector & x)
        {
            // x = values(4, 0.0);
            // x = values(4, 10.0);
            x = values(4, 25.0);
        }

    };

}

#endif //UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP
