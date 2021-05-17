#ifndef UTOPIA_SIMPLE_NEWTON_HPP
#define UTOPIA_SIMPLE_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"

#include "utopia_ConjugateGradient.hpp"
#include "utopia_LinearSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class SimpleNewton final : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;

        SimpleNewton(
            const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>())
            : linear_solver_(linear_solver) {}

        void read(Input &in) override {
            if (linear_solver_) {
                in.get("linear_solver", *linear_solver_);
            } else {
                assert(false);
            }

            in.get("max_it", max_it_);
        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) {
            Matrix H;
            Vector g, inc;

            inc.zeros(layout(x));

            bool converged = false;
            for (int i = 0; i < max_it_; ++i) {
                fun.update(x);
                fun.gradient(x, g);

                Scalar norm_g = norm2(g);

                if (verbose_) utopia::out() << i << ") norm_g: " << norm_g << '\n';

                if (norm_g < 1e-8) {
                    converged = true;
                    break;
                }

                fun.hessian(x, H);
                linear_solver_->solve(H, g, inc);
                x -= inc;

                Scalar norm_inc = norm2(inc);

                if (verbose_) utopia::out() << i << ") norm_inc: " << norm_inc << '\n';

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
        int max_it_{20};
    };

}  // namespace utopia

#endif  // UTOPIA_SIMPLE_NEWTON_HPP
