#ifndef UTOPIA_SMALL_SINGULAR_EXAMPLE_HPP
#define UTOPIA_SMALL_SINGULAR_EXAMPLE_HPP

#include <cassert>
#include <cmath>
#include <vector>
#include "utopia_Function.hpp"

namespace utopia
{
    template<class Matrix, class Vector>
    class SmallSingularExample : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        SmallSingularExample() = default;

        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            assert(x.comm().size() == 1);
            assert(x.size() == 2);
            Vector g;
            gradient(x, g);
            result = 0.5 * norm2(g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(x.comm().size() == 1);

            assert(x.size() == 2);
            g.zeros(layout(x));

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
            assert(x.comm().size() == 1);
            assert(x.size() == 2);

            H.dense(serial_layout(2, 2));

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
}

#endif //UTOPIA_SMALL_SINGULAR_EXAMPLE_HPP
