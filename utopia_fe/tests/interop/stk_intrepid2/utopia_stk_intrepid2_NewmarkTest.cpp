#include "utopia_fe_base.hpp"

#if defined(UTOPIA_ENABLE_INTREPID2) && defined(UTOPIA_ENABLE_TINY_EXPR)

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia.hpp"

#include "utopia_NLSolve.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_VelocityNewmarkIntegrator.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

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
    using Newmark_t = utopia::NewmarkIntegrator<FS_t>;
    using VelocityNewmark_t = utopia::VelocityNewmarkIntegrator<FS_t>;

    std::shared_ptr<FS_t> make_space(const bool with_time_dependent_bc) {
        auto space = std::make_shared<FS_t>();

        auto mesh = param_list(param("type", "cube"), param("nx", 2), param("ny", 2), param("nz", 2));

        if (with_time_dependent_bc) {
            auto params = param_list(
                param("n_var", 1),
                param("mesh", std::move(mesh)),
                param("boundary_conditions",
                      param_list(param("bc0",
                                       param_list(param("type", "TimeDependentCondition"),
                                                  param("name", "left"),
                                                  param("var", 0),
                                                  param("value", "0.1*sin(2*t)"))),
                                 param("bc1", param_list(param("name", "right"), param("var", 0), param("value", 0.))))));
            space->read(params);
        } else {
            auto params = param_list(param("n_var", 1), param("mesh", std::move(mesh)));
            space->read(params);
        }

        return space;
    }

    // Every dof is an independent oscillator x'' + 4 x = cos(3 t) (lumped mass and lumped "Mass" stiffness).
    // x(0) = v(0) = 0, hence a(0) = 1 != 0 (the initial state is NOT in equilibrium).
    template <class Integrator>
    std::shared_ptr<Integrator> make_oscillator(const std::shared_ptr<FS_t> &space, const Scalar_t dt, const int n_steps) {
        auto fun = std::make_shared<Integrator>(space);
        auto params = param_list(
            param("assembly",
                  param_list(param("material", param_list(param("type", "Mass"), param("density", 4.), param("lumped", true))),
                             param("forcing_functions",
                                   param_list(param("f0",
                                                    param_list(param("type", "IncrementalForcingFunction"),
                                                               param("is_symbolic", true),
                                                               param("value", "cos(3*t)"),
                                                               param("component", 0))))))),
            param("time", param_list(param("delta", dt), param("n_steps", n_steps))),
            param("output_path", "newmark_test.e"));
        fun->read(params);
        return fun;
    }

    // Diffusion-like wave problem x'' - lap(x) = 0 with time-dependent Dirichlet data on the left face.
    template <class Integrator>
    std::shared_ptr<Integrator> make_wave(const std::shared_ptr<FS_t> &space, const Scalar_t dt, const int n_steps) {
        auto fun = std::make_shared<Integrator>(space);
        auto params =
            param_list(param("assembly",
                             param_list(param("material", param_list(param("type", "LaplaceOperator"), param("coeff", 1.))))),
                       param("time", param_list(param("delta", dt), param("n_steps", n_steps))),
                       param("output_path", "newmark_test.e"));
        fun->read(params);
        return fun;
    }

    template <class Integrator>
    void solve(const std::shared_ptr<Integrator> &fun) {
        NLSolve<FS_t> nl;
        nl.verbose(false);
        nl.init(fun);
        nl.set_solver(make_exact_newton<Matrix_t, Vector_t>());
        utopia_test_assert(nl.solve());
    }

    Scalar_t oscillator_exact(const Scalar_t t) { return (std::cos(3 * t) - std::cos(2 * t)) / (4. - 9.); }

    template <class Integrator>
    Scalar_t oscillator_error(const int n_steps) {
        const Scalar_t T = 1.;
        auto space = make_space(false);
        auto fun = make_oscillator<Integrator>(space, T / n_steps, n_steps);
        solve(fun);

        utopia_test_assert(approxeq(fun->time()->get(), T, 1e-12));

        Vector_t err = fun->solution();
        err -= Vector_t(layout(err), oscillator_exact(T));
        return norm_infty(err);
    }

    template <class Integrator>
    void check_derivatives(const std::string &name) {
        const Scalar_t dt = 0.1;
        auto space = make_space(true);
        auto fun = make_wave<Integrator>(space, dt, 10);

        Vector_t x;
        space->create_vector(x);
        x.set(0.);
        utopia_test_assert(fun->setup_IVP(x));

        // One artificial step, so that v_n and a_n are non-zero
        Vector_t u = x;
        fill_pattern(u, 0.3, 0.1);
        utopia_test_assert(fun->update_IVP(u));
        utopia_test_assert(norm_infty(fun->velocity()) > 0.);
        utopia_test_assert(norm_infty(fun->acceleration()) > 0.);

        Vector_t p = x, d = x;
        fill_pattern(p, 0.2, 0.5);
        fill_pattern(d, 0.05, 1.3);

        Vector_t g_sep, g_comb, g_p, g_pd;
        Matrix_t H_sep, H_comb;

        utopia_test_assert(fun->gradient(p, g_sep));
        utopia_test_assert(fun->hessian(p, H_sep));
        utopia_test_assert(fun->hessian_and_gradient(p, H_comb, g_comb));

        const Scalar_t g_norm = norm_infty(g_sep);
        const Scalar_t g_diff = norm_infty(Vector_t(g_sep - g_comb));
        const Scalar_t H_diff = norm2(Matrix_t(H_sep - H_comb));

        // The problem is linear, so the Jacobian must reproduce the gradient difference exactly
        Vector_t p_plus_d = p + d;
        utopia_test_assert(fun->gradient(p_plus_d, g_pd));
        Vector_t Hd = H_sep * d;
        space->apply_zero_constraints(Hd);
        Vector_t dg = g_pd - g_sep;
        const Scalar_t dg_norm = norm_infty(dg);
        const Scalar_t fd_abs = norm_infty(Vector_t(dg - Hd));
        const Scalar_t fd_diff = fd_abs / dg_norm;

        if (!space->comm().rank()) {
            utopia::out() << name << ": |g| = " << g_norm << ", |g - g(hessian_and_gradient)| = " << g_diff
                          << ", |H - H(hessian_and_gradient)| = " << H_diff
                          << ", |g(p+d) - g(p) - H d| / |g(p+d) - g(p)| = " << fd_diff << "\n";
        }

        utopia_test_assert(g_diff <= 1e-10 * g_norm);
        utopia_test_assert(H_diff <= 1e-10 * norm2(H_sep));
        utopia_test_assert(fd_diff <= 1e-10);
    }

}  // namespace

// Second-order convergence of both forms (detects the load/BC time lag and a(0) = 0)
void stk_newmark_convergence() {
    const int n_levels = 3;
    const int n_steps_coarse = 16;

    Scalar_t err_x[n_levels], err_v[n_levels];

    for (int l = 0; l < n_levels; ++l) {
        const int n_steps = n_steps_coarse << l;
        err_x[l] = oscillator_error<Newmark_t>(n_steps);
        err_v[l] = oscillator_error<VelocityNewmark_t>(n_steps);
    }

    for (int l = 0; l < n_levels; ++l) {
        const Scalar_t rate_x = l ? std::log2(err_x[l - 1] / err_x[l]) : 0.;
        const Scalar_t rate_v = l ? std::log2(err_v[l - 1] / err_v[l]) : 0.;

        if (!mpi_world_rank()) {
            utopia::out() << "stk_newmark_convergence: n_steps = " << (n_steps_coarse << l)
                          << ", error (displacement form) = " << err_x[l] << ", rate = " << rate_x
                          << ", error (velocity form) = " << err_v[l] << ", rate = " << rate_v << "\n";
        }

        // Both forms are algebraically the same scheme
        utopia_test_assert(std::abs(err_x[l] - err_v[l]) <= 1e-8);

        if (l) {
            utopia_test_assert(rate_x >= 1.9);
            utopia_test_assert(rate_v >= 1.9);
        }
    }
}

UTOPIA_REGISTER_TEST_FUNCTION(stk_newmark_convergence);

// Time-dependent Dirichlet data must be imposed at t_{n+1}, also by the velocity form
void stk_newmark_dirichlet() {
    const Scalar_t dt = 1. / 16;
    const int n_steps = 8;

    auto space_x = make_space(true);
    auto fun_x = make_wave<Newmark_t>(space_x, dt, n_steps);
    solve(fun_x);

    auto space_v = make_space(true);
    auto fun_v = make_wave<VelocityNewmark_t>(space_v, dt, n_steps);
    solve(fun_v);

    const Scalar_t t_end = dt * n_steps;
    const Scalar_t bc_x = dirichlet_violation(*space_x, fun_x->solution());
    const Scalar_t bc_v = dirichlet_violation(*space_v, fun_v->solution());
    const Scalar_t diff = norm_infty(Vector_t(fun_x->solution() - fun_v->solution()));
    const Scalar_t g_end = 0.1 * std::sin(2 * t_end);

    if (!mpi_world_rank()) {
        utopia::out() << "stk_newmark_dirichlet: t = " << fun_x->time()->get() << ", g(t) = " << g_end
                      << ", BC violation (displacement form) = " << bc_x
                      << ", BC violation (velocity form) = " << bc_v << ", |x_disp - x_vel| = " << diff << "\n";
    }

    utopia_test_assert(approxeq(fun_x->time()->get(), t_end, 1e-12));
    utopia_test_assert(approxeq(fun_v->time()->get(), t_end, 1e-12));
    utopia_test_assert(bc_x <= 1e-12);
    utopia_test_assert(bc_v <= 1e-12);
    utopia_test_assert(diff <= 1e-8);
}

UTOPIA_REGISTER_TEST_FUNCTION(stk_newmark_dirichlet);

// gradient / hessian / hessian_and_gradient must describe the same residual and Jacobian
void stk_newmark_derivatives() {
    check_derivatives<Newmark_t>("stk_newmark_derivatives(NewmarkIntegrator)");
    check_derivatives<VelocityNewmark_t>("stk_newmark_derivatives(VelocityNewmarkIntegrator)");
}

UTOPIA_REGISTER_TEST_FUNCTION(stk_newmark_derivatives);

// The default (lumped) mass matrix of a vector-valued space is diagonal, symmetric, positive and integrates
// the density over the domain once per component
void stk_lumped_mass() {
    for (int n_var = 1; n_var <= 3; ++n_var) {
        auto space = std::make_shared<FS_t>();
        auto params =
            param_list(param("n_var", n_var),
                       param("mesh", param_list(param("type", "cube"), param("nx", 2), param("ny", 2), param("nz", 2))));
        space->read(params);

        FEModelFunction<FS_t> fun(space);
        auto assembly = param_list(
            param("assembly",
                  param_list(param("material", param_list(param("type", "LaplaceOperator"), param("coeff", 1.))))));
        fun.read(assembly);
        utopia_test_assert(fun.assemble_mass_matrix());

        const Matrix_t &M = *fun.mass_matrix();
        Vector_t d = diag(M);
        Matrix_t D = diag(d);
        Matrix_t Mt = transpose(M);

        const Scalar_t off_diagonal = norm_infty(Matrix_t(M - D));
        const Scalar_t asymmetry = norm_infty(Matrix_t(M - Mt));
        const Scalar_t total = sum(M);
        const Scalar_t min_diag = min(d);

        if (!space->comm().rank()) {
            utopia::out() << "stk_lumped_mass: n_var = " << n_var << ", |M - diag(M)| = " << off_diagonal
                          << ", |M - M^T| = " << asymmetry << ", sum(M) = " << total << ", min(diag(M)) = " << min_diag
                          << "\n";
        }

        utopia_test_assert(off_diagonal == 0.);
        utopia_test_assert(asymmetry == 0.);
        utopia_test_assert(approxeq(total, Scalar_t(n_var), 1e-12));
        utopia_test_assert(min_diag > 0.);
    }
}

UTOPIA_REGISTER_TEST_FUNCTION(stk_lumped_mass);

#endif  // UTOPIA_ENABLE_INTREPID2 && UTOPIA_ENABLE_TINY_EXPR
