#ifndef UTOPIA_SOLVER_RASTRIGIN_TESTFUNCTION_HPP
#define UTOPIA_SOLVER_RASTRIGIN_TESTFUNCTION_HPP

#include <cassert>
#include <cmath>
#include <vector>

#include "utopia_Function.hpp"

namespace utopia {

    /**
     * @brief     Rastrigin's test function has a lot of regularly distributed local minima, but just one global min. \n
     *            GLOBAL min: \f$ f(x^*) = 0 \text{at} x^* = (0, 0, 0, ..., 0). \n
     *            Function has following structure: \n
     *            \f$ f(x) = 10 * d + \sum_{i = 1}^{d} [ x_i^2 - 10 * cos(2 pi x_i)]  \f$
     *
     *
     */
    template <class Matrix, class Vector>
    class Rastrigin : public Function<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        Rastrigin() : pi(3.141592) {
            help_1 = make_unique<Vector>();
            help_2 = make_unique<Vector>();
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            // *help_1 = pow2(point) - 10. * cos((2. * pi) * point);

            *help_1 = pow2(point);
            *help_2 = (2. * pi) * point;
            *help_2 = cos(*help_2);
            *help_1 = *help_1 - 10. * *help_2;

            result = 10. * point.size() + sum(*help_1);
            return true;
        }

        bool gradient(const Vector &point, Vector &gradient) const override {
            using std::sin;

            if (empty(gradient)) {
                gradient.zeros(layout(point));
            }

            const Read<Vector> read(point);
            const Write<Vector> write(gradient);
            const Range rr = range(point);
            const SizeType r_begin = rr.begin();
            const SizeType r_end = rr.end();

            for (SizeType i = r_begin; i != r_end; i++) {
                const auto x = point.get(i);
                gradient.set(i, 2. * x + 20. * pi * sin(2. * pi * x));
            }

            return true;
        }

        bool hessian(const Vector &point, Matrix &hessian) const override {
            using std::cos;

            // const auto n = point.size();

            if (empty(hessian)) {
                hessian.sparse(square_matrix_layout(layout(point)), 1, 0.0);
            }

            const Read<Vector> read(point);
            const Write<Matrix> write(hessian);

            Range rr = range(point);
            const SizeType r_begin = rr.begin();
            const SizeType r_end = rr.end();

            for (SizeType i = r_begin; i != r_end; i++) {
                const auto x = point.get(i);
                hessian.set(i, i, 40 * pi * pi * cos(2. * pi * x) + 2);
            }

            return true;
        }

    private:
        Scalar pi;
        std::unique_ptr<Vector> help_1;
        std::unique_ptr<Vector> help_2;
    };
}  // namespace utopia

#endif  // UTOPIA_SOLVER_RASTRIGIN_TESTFUNCTION_HPP
