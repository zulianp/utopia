#include "utopia_fe_base.hpp"

#if defined(UTOPIA_ENABLE_INTREPID2) && defined(UTOPIA_ENABLE_TINY_EXPR)

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia.hpp"

#include "utopia_NLSolve.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"

#include "utopia_ObstacleVelocityNewmark.hpp"

#include "utopia_NewmarkTestUtils.hpp"

#include <cmath>
#include <memory>

using namespace utopia;
using namespace utopia::newmark_test;

namespace {

    using FS_t = utopia::stk::FunctionSpace;
    using Matrix_t = Traits<FS_t>::Matrix;
    using Vector_t = Traits<FS_t>::Vector;
    using Scalar_t = Traits<FS_t>::Scalar;
    using ObstacleVelocityNewmark_t = utopia::ObstacleVelocityNewmark<FS_t>;

    // Unit cube, the left face is moved in x-direction by g(t) = 0.05 sin(2 t)
    std::shared_ptr<FS_t> make_space() {
        auto space = std::make_shared<FS_t>();

        auto params = param_list(
            param("n_var", 3),
            param("mesh", param_list(param("type", "cube"), param("nx", 2), param("ny", 2), param("nz", 2))),
            param("boundary_conditions",
                  param_list(param("bc0",
                                   param_list(param("type", "TimeDependentCondition"),
                                              param("name", "left"),
                                              param("var", 0),
                                              param("value", "0.05*sin(2*t)"))),
                             param("bc1", param_list(param("name", "left"), param("var", 1), param("value", 0.))),
                             param("bc2", param_list(param("name", "left"), param("var", 2), param("value", 0.))))));

        space->read(params);
        return space;
    }

    // Plane obstacle at x = 1 + gap in front of the right face. The barrier is active for gaps below
    // its thickness (0.01 by default)
    std::shared_ptr<ObstacleVelocityNewmark_t> make_obstacle_problem(const std::shared_ptr<FS_t> &space,
                                                                     const Scalar_t gap,
                                                                     const Scalar_t dt,
                                                                     const int n_steps) {
        auto fun = std::make_shared<ObstacleVelocityNewmark_t>(std::make_shared<FEModelFunction<FS_t>>(space));

        auto params = param_list(
            param("assembly",
                  param_list(param("material",
                                   param_list(param("type", "LinearElasticity"), param("mu", 1.), param("lambda", 1.))))),
            param("time", param_list(param("delta", dt), param("n_steps", n_steps))),
            param("obstacle",
                  param_list(param("type", "analytic"),
                             param("fun",
                                   param_list(param("plane",
                                                    param_list(param("type", "plane"),
                                                               param("side", "right"),
                                                               param("d", -(1 + gap)),
                                                               param("normal_x", -1.),
                                                               param("normal_y", 0.),
                                                               param("normal_z", 0.))))))),
            param("function_type", "BoundedLogBarrier"),
            param("barrier_parameter", 1e-3),
            param("min_barrier_parameter", 1e-3),
            param("output_path", "newmark_obstacle_test.e"));

        fun->read(params);
        return fun;
    }

    // Problem after setup_IVP and one artificial step, so that v_n and a_n are non-zero
    std::shared_ptr<ObstacleVelocityNewmark_t> make_stepped_obstacle_problem(const std::shared_ptr<FS_t> &space,
                                                                             const Scalar_t gap,
                                                                             const Scalar_t dt,
                                                                             const int n_steps) {
        auto fun = make_obstacle_problem(space, gap, dt, n_steps);

        Vector_t x;
        space->create_vector(x);
        x.set(0.);
        utopia_test_assert(fun->setup_IVP(x));

        Vector_t u = x;
        fill_pattern(u, 0.01, 0.1);
        utopia_test_assert(fun->update_IVP(u));
        return fun;
    }

}  // namespace

// ObstacleVelocityNewmark: gradient, hessian and hessian_and_gradient describe the same residual and Jacobian
// (including the barrier), and time-dependent Dirichlet data is imposed at t_{n+1}
void moonolith_stk_newmark_obstacle() {
    const Scalar_t dt = 1. / 32;
    const int n_steps = 4;

    {
        // The obstacle is closer than the barrier thickness, so the barrier contributes
        auto space = make_space();
        auto fun = make_stepped_obstacle_problem(space, 0.005, dt, n_steps);

        Vector_t p, d;
        space->create_vector(p);
        space->create_vector(d);
        fill_pattern(p, 0.01, 0.5);
        fill_pattern(d, 0.01, 1.3);

        // As in Newton
        utopia_test_assert(fun->update(p));

        Vector_t g_sep, g_comb, g_plus, g_minus;
        Matrix_t H_sep, H_comb;

        utopia_test_assert(fun->gradient(p, g_sep));
        utopia_test_assert(fun->hessian(p, H_sep));
        utopia_test_assert(fun->hessian_and_gradient(p, H_comb, g_comb));

        const Scalar_t g_norm = norm_infty(g_sep);
        const Scalar_t g_diff = norm_infty(Vector_t(g_sep - g_comb));
        const Scalar_t H_diff = norm2(Matrix_t(H_sep - H_comb));

        // Same state with a far obstacle (inactive barrier)
        auto space_far = make_space();
        auto fun_far = make_stepped_obstacle_problem(space_far, 0.5, dt, n_steps);
        utopia_test_assert(fun_far->update(p));
        Vector_t g_far;
        utopia_test_assert(fun_far->gradient(p, g_far));
        const Scalar_t barrier_contribution = norm_infty(Vector_t(g_sep - g_far));

        // Central finite difference of the gradient (the barrier makes the residual nonlinear)
        const Scalar_t eps = 1e-6;
        Vector_t p_plus = p + eps * d;
        Vector_t p_minus = p - eps * d;
        utopia_test_assert(fun->gradient(p_plus, g_plus));
        utopia_test_assert(fun->gradient(p_minus, g_minus));

        Vector_t dg = (1. / (2 * eps)) * (g_plus - g_minus);
        Vector_t Hd = H_sep * d;
        space->apply_zero_constraints(Hd);

        const Scalar_t dg_norm = norm_infty(dg);
        const Scalar_t fd_abs = norm_infty(Vector_t(dg - Hd));
        const Scalar_t fd_diff = fd_abs / dg_norm;

        if (!space->comm().rank()) {
            utopia::out() << "moonolith_stk_newmark_obstacle: |g| = " << g_norm
                          << ", |barrier contribution to g| = " << barrier_contribution
                          << ", |g - g(hessian_and_gradient)| = " << g_diff
                          << ", |H - H(hessian_and_gradient)| = " << H_diff << ", |FD(g) d - H d| / |FD(g) d| = " << fd_diff
                          << "\n";
        }

        utopia_test_assert(barrier_contribution > 1e-6 * g_norm);
        utopia_test_assert(g_diff <= 1e-10 * g_norm);
        utopia_test_assert(H_diff <= 1e-10 * norm2(H_sep));
        utopia_test_assert(fd_diff <= 1e-5);
    }

    {
        auto space = make_space();
        auto fun = make_obstacle_problem(space, 0.05, dt, n_steps);

        NLSolve<FS_t> nl;
        nl.verbose(false);
        nl.init(fun);
        nl.set_solver(make_exact_newton<Matrix_t, Vector_t>());
        utopia_test_assert(nl.solve());

        const Scalar_t t_end = dt * n_steps;
        const Scalar_t bc = dirichlet_violation(*space, fun->solution());

        if (!space->comm().rank()) {
            utopia::out() << "moonolith_stk_newmark_obstacle: t = " << fun->time()->get()
                          << ", g(t) = " << 0.05 * std::sin(2 * t_end) << ", BC violation = " << bc << "\n";
        }

        utopia_test_assert(approxeq(fun->time()->get(), t_end, 1e-12));
        utopia_test_assert(bc <= 1e-12);
    }
}

UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_obstacle);

#endif  // UTOPIA_ENABLE_INTREPID2 && UTOPIA_ENABLE_TINY_EXPR
