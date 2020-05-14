#ifndef UTOPIA_QUASI_NEWTON_BOUND_HPP
#define UTOPIA_QUASI_NEWTON_BOUND_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_QuasiNewtonBase.hpp"
#include "utopia_Layout.hpp"

#include <iomanip>
#include <limits>

namespace utopia
{
    template<class Vector>
    class QuasiNewtonBound : public QuasiNewtonBase<Vector>,
                             public VariableBoundSolverInterface<Vector>

    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        using LSStrategy           = utopia::LSStrategy<Vector>;
        using HessianApproximation = utopia::HessianApproximation<Vector>;
        using QPSolver             = utopia::MatrixFreeQPSolver<Vector>;

    public:

        QuasiNewtonBound(
            const std::shared_ptr <HessianApproximation> &hessian_approx,
            const std::shared_ptr <QPSolver> &linear_solver
        ):
            QuasiNewtonBase<Vector>(hessian_approx, linear_solver),
            initialized_(false)
        {}

        bool solve(FunctionBase<Vector> &fun, Vector &x) override
        {
            using namespace utopia;

            Scalar g_norm, r_norm, g0_norm, s_norm = 1;
            SizeType it = 0;

            Scalar alpha = 1.0;
            bool converged = false;

            const auto x_layout = utopia::layout(x);

            this->fill_empty_bounds(x_layout);
            this->make_iterate_feasible(x);

            if(!initialized_ || !x_layout.same(layout_))
            {
                init_memory(x_layout);
            }

            fun.gradient(x, g);
            g0_norm = this->criticality_measure_infty(x, g);

            QuasiNewtonBase<Vector>::init_memory(x, g);

            if(this->verbose_) {
                this->init_solver("QUASI NEWTON BOUND", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g0_norm, s_norm});
            }

            it++;

            this->initialize_approximation(x, g);
            auto multiplication_action = this->hessian_approx_strategy_->build_apply_H();

            UTOPIA_NO_ALLOC_BEGIN("Quasi_NewtonBound");

            while(!converged)
            {
                if(QPSolver * qp_solver = dynamic_cast<QPSolver*>(this->linear_solver().get()))
                {
                    auto box = this->build_correction_constraints(x);
                    qp_solver->set_box_constraints(box);
                    s.set(0.0);
                    g_minus = -1.0 * g;
                    qp_solver->solve(*multiplication_action, g_minus, s);
                }
                else
                {
                    utopia_error("utopia::QuasiNewtonBound: MF solver which is not QPSolver is not suported at the moment... \n");
                }

                // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton Bound :3");
                alpha = this->get_alpha(fun, g, x, s);
                s *= alpha;

                x = x + s;
                this->make_iterate_feasible(x);

                y = g;
                // UTOPIA_NO_ALLOC_END();


                fun.gradient(x, g);

                // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton Bound :4");
                // norms needed for convergence check
                g_norm = this->criticality_measure_infty(x, g);
                s_norm = norm_infty(s);
                r_norm = g_norm/g0_norm;

                // diff between fresh and old grad...
                y = g - y;
                this->update(s, y, x, g);
                // UTOPIA_NO_ALLOC_END();

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm,  s_norm, alpha});

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);

                it++;
            }

            UTOPIA_NO_ALLOC_END();

            this->print_statistics(it);
            return true;
        }

    private:
        void init_memory(const Layout &layout) override {
            s.zeros(layout);
            g.zeros(layout);
            y.zeros(layout);
            g_minus.zeros(layout);

            VariableBoundSolverInterface<Vector>::init_memory(layout);

            initialized_ = true;
            layout_ = layout;
        }

        Vector g, s, y, g_minus;
        bool initialized_;
        Layout layout_;
    };
}

#endif //UTOPIA_QUASI_NEWTON_BOUND_HPP
