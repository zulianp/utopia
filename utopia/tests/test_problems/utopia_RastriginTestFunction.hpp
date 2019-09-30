#ifndef UTOPIA_SOLVER_RASTRIGIN_TESTFUNCTION_HPP
#define UTOPIA_SOLVER_RASTRIGIN_TESTFUNCTION_HPP

#include <vector>
#include <cassert>
#include <cmath>

#include "utopia_Function.hpp"


namespace utopia
{

    /**
     * @brief     Rastrigin's test function has a lot of regularly distributed local minima, but just one global min. \n
     *            GLOBAL min: \f$ f(x^*) = 0 \text{at} x^* = (0, 0, 0, ..., 0). \n
     *            Function has following structure: \n
     *            \f$ f(x) = 10 * d + \sum_{i = 1}^{d} [ x_i^2 - 10 * cos(2 pi x_i)]  \f$
     *
     *
     */
    template<class Matrix, class Vector>
    class Rastrigin : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        Rastrigin() : pi(3.141592) {}

        bool value(const Vector &point, typename Vector::Scalar &result) const override
        {
            result = 10. * point.size() + sum(pow2(point) - 10. * cos((2. * pi) * point));
            return true;
        }

        bool gradient(const Vector &point, Vector &gradient) const override
        {
            using std::sin;

            if(empty(gradient)) {
                gradient = zeros(point.size());
            }

            const Read<Vector> read(point);
            const Write<Vector> write(gradient);
            const Range rr = range(point);

            for(SizeType i = rr.begin(); i != rr.end(); i++) {
                const auto x = point.get(i);
                gradient.set(i, 2. * x + 20. * pi * sin(2. * pi * x));
            }

            return true;
        }

        bool hessian(const Vector &point, Matrix &hessian) const override
        {
            using std::cos;

            const auto n = point.size();

            if(empty(hessian)) {
                hessian = sparse(n, n, 1);
            }

            const Read<Vector> read(point);
            const Write<Matrix> write(hessian);

            Range rr = range(point);

            for(SizeType i = rr.begin(); i != rr.end(); i++) {
                const auto x = point.get(i);
                hessian.set(i, i,  40 * pi * pi * cos(2. * pi * x) + 2 );
            }

            return true;
        }

    private:
        Scalar pi;
    };
}

#endif //UTOPIA_SOLVER_RASTRIGIN_TESTFUNCTION_HPP
