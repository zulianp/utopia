#ifndef UTOPIA_SIMPLE_NEWTON_HPP
#define UTOPIA_SIMPLE_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"

#include "utopia_ConjugateGradient.hpp"
#include "utopia_LinearSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class SimpleNewton final {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;

        SimpleNewton(
            const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>())
            : linear_solver_(linear_solver) {}

        bool solve(Function<Matrix, Vector> &fun, Vector &x) {
            Matrix H;
            Vector g, inc;

            inc.zeros(layout(x));

            bool converged = false;
            for (int i = 0; i < 20; ++i) {
                // function_->hessian_and_gradient(x, H, g);
                fun.update(x);
                fun.hessian_and_gradient(x, H, g);

                Scalar norm_g = norm2(g);

                if (verbose_) utopia::out() << "norm_g: " << norm_g << '\n';

                if (norm_g < 1e-8) {
                    converged = true;
                    break;
                }

                linear_solver_->solve(H, g, inc);
                x -= inc;
                Scalar norm_inc = norm2(inc);

                if (verbose_) utopia::out() << "norm_inc: " << norm_inc << '\n';

                if (norm_inc < 1e-8) {
                    converged = true;
                    break;
                }
                inc.set(0.0);
            }

            return true;
        }

        inline void verbose(const bool verbose) { verbose_ = verbose; }

    private:
        std::shared_ptr<LinearSolver> linear_solver_;
        bool verbose_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_SIMPLE_NEWTON_HPP
