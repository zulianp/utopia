#ifndef UTOPIA_NEWMARK_TEST_UTILS_HPP
#define UTOPIA_NEWMARK_TEST_UTILS_HPP

#include "utopia.hpp"

#include <cmath>
#include <memory>

namespace utopia {
    namespace newmark_test {

        // Newton with an exact (direct) linear solve, so that the discrete equations are satisfied
        // to round-off and only the time discretisation error remains.
        template <class Matrix, class Vector>
        std::shared_ptr<Newton<Matrix, Vector>> make_exact_newton() {
            auto newton = std::make_shared<Newton<Matrix, Vector>>(std::make_shared<Factorization<Matrix, Vector>>());
            newton->atol(1e-11);
            newton->rtol(1e-14);
            newton->stol(1e-14);
            newton->max_it(10);
            newton->verbose(false);
            return newton;
        }

        // Fills v with a smooth, non-constant pattern
        template <class Vector>
        void fill_pattern(Vector &v, const typename Traits<Vector>::Scalar scale, const typename Traits<Vector>::Scalar phase) {
            using SizeType = typename Traits<Vector>::SizeType;

            auto r = range(v);
            auto view = local_view_device(v);
            const SizeType offset = r.begin();
            parallel_for(
                local_range_device(v), UTOPIA_LAMBDA(const SizeType i) {
                    view.set(i, scale * std::sin(phase + 0.7 * (offset + i)));
                });
        }

        // Distance of the constrained dofs of x from the Dirichlet data at the current time of the space
        template <class FunctionSpace, class Vector>
        typename Traits<Vector>::Scalar dirichlet_violation(const FunctionSpace &space, const Vector &x) {
            Vector bc = x;
            space.apply_constraints(bc);
            bc -= x;
            return norm_infty(bc);
        }

    }  // namespace newmark_test
}  // namespace utopia

#endif  // UTOPIA_NEWMARK_TEST_UTILS_HPP
