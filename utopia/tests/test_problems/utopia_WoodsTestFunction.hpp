/*
* @Author: Alena Kopanicakova
* @Date:   2016-08-01
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2016-08-01
*/
#ifndef UTOPIA_SOLVER_WOODS_TESTFUNCTION_HPP
#define UTOPIA_SOLVER_WOODS_TESTFUNCTION_HPP

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia 
{

    /**
     * @brief       Wood's test function widely used for testing nonlinear optimization problems. \n
     *              It is a fourth-degree polynomial which is reasonably well-behaved near the minimum,
     *               but in order to get there one must cross a rather flat, four-dimensional ‘plateau’
     *               which often causes minimization algorithm to get ‘stuck’ far from the minimum. As
     *               such it is a particularly good test of convergence criteria and simulates quite well a
     *               feature of many physical problems in many variables where no good starting approximation
     *               is known. \n
     *               Exact solution is at F(1, 1, 1, 1) = 0. \n
     *               Good intial guess for testing: F(-3, -1, -3, -1) = 19192. \n
     */
    template<class Matrix, class Vector>
    class Woods : public Function<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        Woods() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == 4);

            {
                const Read<Vector> read(point);

                const Scalar w = point.get(0);
                const Scalar x = point.get(1);
                const Scalar y = point.get(2);
                const Scalar z = point.get(3);

                result = 100* pow( x - w*w, 2) + pow(w - 1, 2) + 90 * pow(z -y*y, 2) + pow(1 -y, 2) + 10.1 * ( pow(x - 1, 2) + pow(z -1, 2)) + 19.8 * (x-1) * (z-1); 
            }
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override 
        {
            assert(point.size().get(0) == 4);
            result = zeros(4);

            {
                const Read<Vector> read(point);
                const Write<Vector> write(result);

                const Scalar w = point.get(0);
                const Scalar x = point.get(1);
                const Scalar y = point.get(2);
                const Scalar z = point.get(3);

                result.set(0,    2 * w - 400 * w * ( -w * w + x) -2);
                result.set(1,   -200 * w* w + (1101 * x)/5 + (99 * z)/5 - 40);
                result.set(2,   2 * y - 360 * y * ( - y * y + z ) - 2 );
                result.set(3,   -180 * y * y + (99 * x)/5 + (1001 * z)/5 -40);
            }
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override 
        {
            assert(point.size().get(0) == 4);
            result = zeros(4, 4);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar w = point.get(0);
            const Scalar x = point.get(1);
            const Scalar y = point.get(2);
            const Scalar z = point.get(3);

            result.set(0, 0, 1200 * w * w - 400 * x +2);
            result.set(0, 1, -400 * w);
            result.set(0, 2, 0);
            result.set(0, 3, 0);

            result.set(1, 0, -400 * w);
            result.set(1, 1, 1101/5);
            result.set(1, 2, 0);
            result.set(1, 3, 99/5);

            result.set(2, 0, 0);
            result.set(2, 1, 0);
            result.set(2, 2, 1080 * y * y - 360 * z + 2 );
            result.set(2, 3, - 360 * y);

            result.set(3, 0, 0);
            result.set(3, 1, 99/5);
            result.set(3, 2, -360 * y);
            result.set(3, 3, 1001/5);
            
            return true;
        }
    };
}

#endif //UTOPIA_SOLVER_WOODS_TESTFUNCTION_HPP
