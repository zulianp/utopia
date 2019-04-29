#ifndef UTOPIA_SOLVER_QUASI_NEWTON_HPP
#define UTOPIA_SOLVER_QUASI_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_QuasiNewtonBase.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Vector>
    class QuasiNewton : public QuasiNewtonBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                               Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                            SizeType;

        typedef utopia::HessianApproximation<Vector>                HessianApproximation;
        typedef utopia::MatrixFreeLinearSolver<Vector>              LinSolver;

    public:

        QuasiNewton(const std::shared_ptr <HessianApproximation> &hessian_approx,
                    const std::shared_ptr <LinSolver> &linear_solver):
                    QuasiNewtonBase<Vector>(hessian_approx, linear_solver)
        {

        }

        bool solve(FunctionBase<Vector> &fun, Vector &x) override
        {
            using namespace utopia;

            Vector g, s, y;

            Scalar g_norm, g0_norm, r_norm=1, s_norm=1;
            SizeType it = 0;

            Scalar alpha = 1.0;

            bool converged = false;

            fun.gradient(x, g);
            g0_norm = norm2(g);
            g_norm = g0_norm;

            if(this->verbose_) {
                this->init_solver("QUASI NEWTON", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm});
            }
            it++;

            this->initialize_approximation();

            while(!converged)
            {
                s = local_zeros(local_size(x));
                this->linear_solve(-1.0 * g, s);

                alpha = this->get_alpha(fun, g, x, s);

                s *= alpha;
                x+=s;

                y = g;
                fun.gradient(x, g);

                // norms needed for convergence check
                norms2(g, s, g_norm, s_norm);
                r_norm = g_norm/g0_norm;

                // diff between fresh and old grad...
                y = g - y;
                this->update(s, y);

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, alpha});

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);

                it++;
            }

            this->print_statistics(it);
            return true;
        }


    };

}
#endif //UTOPIA_SOLVER_QUASI_NEWTON_HPP
