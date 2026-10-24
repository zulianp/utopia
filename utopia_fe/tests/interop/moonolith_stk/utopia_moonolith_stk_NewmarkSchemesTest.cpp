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

#include "utopia_ContactFactory.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_ObstacleFEFunction.hpp"
#include "utopia_ObstacleNewmark.hpp"
#include "utopia_ObstacleStabilizedNewmark.hpp"
#include "utopia_ObstacleStabilizedVelocityNewmark.hpp"
#include "utopia_ObstacleVelocityNewmark.hpp"
#include "utopia_VelocityNewmarkIntegrator.hpp"

#include "utopia_NewmarkTestUtils.hpp"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace utopia;
using namespace utopia::newmark_test;

namespace {

    using FS_t = utopia::stk::FunctionSpace;
    using Matrix_t = Traits<FS_t>::Matrix;
    using Vector_t = Traits<FS_t>::Vector;
    using Scalar_t = Traits<FS_t>::Scalar;
    using SizeType_t = Traits<FS_t>::SizeType;
    using Contact_t = utopia::ContactInterface<FS_t>;

    using Newmark_t = utopia::NewmarkIntegrator<FS_t>;
    using VelocityNewmark_t = utopia::VelocityNewmarkIntegrator<FS_t>;
    using ObstacleNewmark_t = utopia::ObstacleNewmark<FS_t>;
    using ObstacleVelocityNewmark_t = utopia::ObstacleVelocityNewmark<FS_t>;
    using StabilizedNewmark_t = utopia::ObstacleStabilizedNewmark<FS_t>;
    using StabilizedVelocityNewmark_t = utopia::ObstacleStabilizedVelocityNewmark<FS_t>;

    // Linear elasticity on the unit cube with a plane obstacle at x = 1 + gap, in front of the right face
    struct Config {
        Scalar_t dt{1. / 32};
        int n_steps{8};
        Scalar_t gap{0.5};
        // Contact enforcement: BoundedLogBarrier (barrier_parameter) or ShiftedPenalty (penalty_parameter)
        std::string function_type{"BoundedLogBarrier"};
        Scalar_t barrier_parameter{1e-3};
        Scalar_t penalty_parameter{1e3};
        // Tilt of the obstacle plane: normal (-1, normal_y, 0) / |.|, with the smallest gap at y = 0
        Scalar_t normal_y{0};
        // Curved obstacle instead of the plane: sphere of radius sphere_radius, at distance gap from the center of the
        // right face (the normals differ at every node)
        Scalar_t sphere_radius{0};
        Scalar_t barrier_thickness{0.01};
        // Lower bound of the barrier parameter, which shrinks during the Newton iterations (barrier_parameter if < 0)
        Scalar_t min_barrier_parameter{-1};
        // Uniform body force in x-direction (towards the obstacle)
        Scalar_t body_force{0.5};
        bool trivial_obstacle{false};
        bool enable_line_search{false};
        // Step damping of the barrier schemes, not set if <= 0 (default of the scheme)
        Scalar_t damping{0.98};
        // The left face is moved in x-direction by g(t) = 0.05 sin(2 t), and fixed in y and z
        bool with_bc{true};
    };

    std::shared_ptr<FS_t> make_space(const bool with_bc) {
        auto space = std::make_shared<FS_t>();

        auto mesh = param_list(param("type", "cube"), param("nx", 2), param("ny", 2), param("nz", 2));

        if (with_bc) {
            auto params = param_list(
                param("n_var", 3),
                param("mesh", std::move(mesh)),
                param("boundary_conditions",
                      param_list(param("bc0",
                                       param_list(param("type", "TimeDependentCondition"),
                                                  param("name", "left"),
                                                  param("var", 0),
                                                  param("value", "0.05*sin(2*t)"))),
                                 param("bc1", param_list(param("name", "left"), param("var", 1), param("value", 0.))),
                                 param("bc2", param_list(param("name", "left"), param("var", 2), param("value", 0.))))));
            space->read(params);
        } else {
            auto params = param_list(param("n_var", 3), param("mesh", std::move(mesh)));
            space->read(params);
        }

        return space;
    }

    // Plane obstacle n . p = d in front of the right face, with n = (-1, normal_y, 0) / |.|. The gap of the right
    // face is gap + normal_y y / |.| (gap at y = 0)
    InputParameters obstacle_params(const Scalar_t gap, const Scalar_t normal_y = 0) {
        const Scalar_t length = std::sqrt(1 + normal_y * normal_y);
        return param_list(param("type", "analytic"),
                          param("fun",
                                param_list(param("plane",
                                                 param_list(param("type", "plane"),
                                                            param("side", "right"),
                                                            param("d", -1. / length - gap),
                                                            param("normal_x", -1. / length),
                                                            param("normal_y", normal_y / length),
                                                            param("normal_z", 0.))))));
    }

    // Sphere of radius R with center (1 + gap + R, 0.5, 0.5). Gap g = |p - c| - R, normal (c - p) / |p - c| (towards
    // the obstacle)
    InputParameters sphere_obstacle_params(const Scalar_t gap, const Scalar_t R) {
        const std::string cx = std::to_string(1 + gap + R);
        const std::string r = "sqrt((x - " + cx + ")^2 + (y - 0.5)^2 + (z - 0.5)^2)";
        return param_list(param("type", "analytic"),
                          param("fun",
                                param_list(param("sphere",
                                                 param_list(param("type", "gap"),
                                                            param("side", "right"),
                                                            param("d", std::to_string(R) + " - " + r),
                                                            param("normal_x", "(" + cx + " - x) / " + r),
                                                            param("normal_y", "(0.5 - y) / " + r),
                                                            param("normal_z", "(0.5 - z) / " + r))))));
    }

    InputParameters obstacle_params(const Config &c) {
        return c.sphere_radius > 0 ? sphere_obstacle_params(c.gap, c.sphere_radius) : obstacle_params(c.gap, c.normal_y);
    }

    InputParameters make_params(const Config &c) {
        auto params = param_list(
            param("assembly",
                  param_list(param("material",
                                   param_list(param("type", "LinearElasticity"), param("mu", 1.), param("lambda", 1.))),
                             param("forcing_functions",
                                   param_list(param("f0",
                                                    param_list(param("type", "IncrementalForcingFunction"),
                                                               param("is_symbolic", true),
                                                               param("value", std::to_string(c.body_force)),
                                                               param("component", 0))))))),
            param("time", param_list(param("delta", c.dt), param("n_steps", c.n_steps))),
            param("obstacle", obstacle_params(c)),
            param("barrier_thickness", c.barrier_thickness),
            param("function_type", c.function_type),
            param("barrier_parameter", c.barrier_parameter),
            param("penalty_parameter", c.penalty_parameter),
            param("min_barrier_parameter",
                  c.min_barrier_parameter >= 0 ? c.min_barrier_parameter : c.barrier_parameter),
            param("output_path", "newmark_schemes_test.e"));

        params.set("trivial_obstacle", c.trivial_obstacle);
        params.set("enable_line_search", c.enable_line_search);
        if (c.damping > 0) {
            params.set("damping", c.damping);
        }
        return params;
    }

    std::shared_ptr<Contact_t> make_obstacle(const Scalar_t gap) {
        auto params = obstacle_params(gap);
        return ContactFactory<FS_t>::new_obstacle(params);
    }

    // ObstacleStabilizedNewmark gets the obstacle from the contact formulation (e.g., ObstacleFEFunction)
    template <class Scheme>
    void attach_obstacle(Scheme &, const Config &) {}

    void attach_obstacle(StabilizedNewmark_t &fun, const Config &c) {
        auto params = obstacle_params(c);
        fun.set_contact(ContactFactory<FS_t>::new_obstacle(params));
    }

    template <class Scheme>
    std::shared_ptr<Scheme> make_scheme(const std::shared_ptr<FS_t> &space, const Config &c) {
        auto fun = std::make_shared<Scheme>(std::make_shared<FEModelFunction<FS_t>>(space));
        auto params = make_params(c);
        fun->read(params);
        attach_obstacle(*fun, c);
        return fun;
    }

    // Initial displacement that is not in equilibrium and satisfies the Dirichlet data at t = 0
    void initial_displacement(const FS_t &space, Vector_t &x0) {
        space.create_vector(x0);
        fill_pattern(x0, 0.002, 0.3);
        space.apply_zero_constraints(x0);
    }

    std::shared_ptr<Newton<Matrix_t, Vector_t>> make_newton() {
        auto newton = make_exact_newton<Matrix_t, Vector_t>();
        // Damped (line-search) steps converge linearly
        newton->max_it(200);
        newton->atol(1e-10);
        return newton;
    }

    template <class Scheme>
    Vector_t run(const Config &c, Scalar_t &bc_violation) {
        auto space = make_space(c.with_bc);
        auto fun = make_scheme<Scheme>(space, c);

        auto x0 = std::make_shared<Vector_t>();
        initial_displacement(*space, *x0);

        NLSolve<FS_t> nl;
        nl.verbose(false);
        nl.init(fun);
        nl.set_solution(std::make_shared<Field<FS_t>>("x", space, x0));
        nl.set_solver(make_newton());
        utopia_test_assert(nl.solve());
        utopia_test_assert(approxeq(fun->time()->get(), c.dt * c.n_steps, 1e-12));

        bc_violation = dirichlet_violation(*space, fun->solution());
        return fun->solution();
    }

    // Vector with value s on the x-components and 0 elsewhere
    void fill_x_component(Vector_t &v, const Scalar_t s) {
        auto r = range(v);
        auto view = local_view_device(v);
        const SizeType_t offset = r.begin();
        parallel_for(
            local_range_device(v), UTOPIA_LAMBDA(const SizeType_t i) { view.set(i, ((offset + i) % 3 == 0) ? s : 0.); });
    }

    // Largest violation max(T(increment) - gap) of the obstacle constraint, with the gap measured in configuration
    // x_ref. The gap is "infinite" for the components that are not in contact
    Scalar_t constraint_violation(FS_t &space,
                                  Contact_t &obstacle,
                                  const Vector_t &x_ref,
                                  const Vector_t &increment) {
        space.displace(x_ref);
        obstacle.assemble(space);
        space.displace(-x_ref);

        Vector_t t_inc;
        obstacle.transform(increment, t_inc);
        t_inc -= obstacle.gap();
        return max(t_inc);
    }

    template <class Scheme>
    void barrier_contribution(const Config &c,
                              const Vector_t &u,
                              const Vector_t &p,
                              Vector_t &g_barrier,
                              Vector_t &g,
                              Matrix_t &H,
                              std::shared_ptr<Scheme> &fun) {
        auto space = make_space(c.with_bc);
        fun = make_scheme<Scheme>(space, c);

        Vector_t x0(layout(u), 0.);
        utopia_test_assert(fun->setup_IVP(x0));
        utopia_test_assert(fun->update_IVP(u));
        utopia_test_assert(fun->gradient(p, g));
        utopia_test_assert(fun->hessian(p, H));

        // Same state without contact term (zero barrier or penalty parameter). The obstacle stays, since the
        // predictor of the stabilized scheme depends on it
        Config c_off = c;
        c_off.barrier_parameter = 0;
        c_off.min_barrier_parameter = 0;
        c_off.penalty_parameter = 0;
        auto space_off = make_space(c.with_bc);
        auto fun_off = make_scheme<Scheme>(space_off, c_off);
        utopia_test_assert(fun_off->setup_IVP(x0));
        utopia_test_assert(fun_off->update_IVP(u));

        Vector_t g_off;
        utopia_test_assert(fun_off->gradient(p, g_off));
        g_barrier = g - g_off;
    }

    // gradient, hessian, hessian_and_gradient and finite differences describe the same residual and Jacobian,
    // with an active barrier (or penalty), and the contact term scales with the configured barrier_parameter
    // (or penalty_parameter). The penalty needs penetration, hence the negative gap. A tilted plane makes the
    // transformation to the normal-tangential frame non-diagonal
    template <class Scheme>
    void check_obstacle_derivatives(const std::string &name,
                                    const Scalar_t u_scale,
                                    const Scalar_t p_scale,
                                    const std::string &function_type = "BoundedLogBarrier",
                                    const Scalar_t normal_y = 0,
                                    const Scalar_t sphere_radius = 0) {
        Config c;
        c.function_type = function_type;
        c.gap = function_type == "ShiftedPenalty" ? -0.002 : 0.005;
        c.normal_y = normal_y;

        if (sphere_radius > 0) {
            // All nodes of the right face in contact: the corners are at distance 0.23 from the sphere (gap 0.005)
            c.sphere_radius = sphere_radius;
            c.barrier_thickness = 0.25;
            if (function_type == "ShiftedPenalty") c.gap = -0.25;
        }

        auto space = make_space(c.with_bc);
        Vector_t u, p, d;
        space->create_vector(u);
        space->create_vector(p);
        space->create_vector(d);

        // Previous step and evaluation point (displacement or velocity, depending on the scheme)
        fill_pattern(u, u_scale, 0.1);
        fill_pattern(p, p_scale, 0.5);
        fill_pattern(d, p_scale, 1.3);

        std::shared_ptr<Scheme> fun, fun_2;
        Vector_t g_barrier, g, g_barrier_2, g_2;
        Matrix_t H, H_2;
        barrier_contribution(c, u, p, g_barrier, g, H, fun);

        Config c_2 = c;
        c_2.barrier_parameter = 2 * c.barrier_parameter;
        c_2.penalty_parameter = 2 * c.penalty_parameter;
        barrier_contribution(c_2, u, p, g_barrier_2, g_2, H_2, fun_2);

        Vector_t g_comb, g_plus, g_minus;
        Matrix_t H_comb;
        utopia_test_assert(fun->hessian_and_gradient(p, H_comb, g_comb));

        const Scalar_t g_norm = norm_infty(g);
        const Scalar_t barrier_norm = norm_infty(g_barrier);
        const Scalar_t g_diff = norm_infty(Vector_t(g - g_comb));
        const Scalar_t H_diff = norm2(Matrix_t(H - H_comb));

        // The configured barrier_parameter (penalty_parameter) is used from the start (no update() is called)
        const Scalar_t scaling_error = norm_infty(Vector_t(g_barrier_2 - 2. * g_barrier)) / barrier_norm;

        // Central finite difference of the gradient (the barrier makes the residual nonlinear)
        const Scalar_t eps = 1e-6;
        Vector_t p_plus = p + eps * d;
        Vector_t p_minus = p - eps * d;
        utopia_test_assert(fun->gradient(p_plus, g_plus));
        utopia_test_assert(fun->gradient(p_minus, g_minus));

        Vector_t dg = (1. / (2 * eps)) * (g_plus - g_minus);
        Vector_t Hd = H * d;
        space->apply_zero_constraints(Hd);
        const Scalar_t dg_norm = norm_infty(dg);
        const Scalar_t fd_diff = norm_infty(Vector_t(dg - Hd)) / dg_norm;

        if (!space->comm().rank()) {
            utopia::out() << name << " (" << function_type
                          << (sphere_radius > 0 ? ", sphere" : (normal_y != 0 ? ", tilted plane" : ""))
                          << "): |g| = " << g_norm << ", |barrier part of g| = " << barrier_norm
                          << ", barrier scaling error = " << scaling_error
                          << ", |g - g(hessian_and_gradient)| = " << g_diff
                          << ", |H - H(hessian_and_gradient)| = " << H_diff
                          << ", |FD(g) d - H d| / |FD(g) d| = " << fd_diff << "\n";
        }

        utopia_test_assert(barrier_norm > 1e-3 * g_norm);
        utopia_test_assert(scaling_error <= 1e-8);
        utopia_test_assert(g_diff <= 1e-10 * g_norm);
        utopia_test_assert(H_diff <= 1e-10 * norm2(H));
        utopia_test_assert(fd_diff <= 1e-5);
    }

    // Step length of the line-search from the configuration x0 = 0.3 gap e_x, for the correction
    // that moves the displacement by gap e_x. The admissible step is 0.7, for both the trivial obstacle
    // (gap of the reference configuration, constraint on x) and the updated one (gap of x0, constraint on x - x0)
    template <class Scheme>
    Scalar_t line_search_alpha(const Config &c, const bool velocity_unknown) {
        auto space = make_space(c.with_bc);
        auto fun = make_scheme<Scheme>(space, c);

        Vector_t x0, g, iterate, correction;
        space->create_vector(x0);
        fill_x_component(x0, 0.3 * c.gap);
        utopia_test_assert(fun->setup_IVP(x0));

        space->create_vector(g);
        space->create_vector(correction);
        g.set(0.);

        if (velocity_unknown) {
            // x(v) = x0 for v = 0 (v_n = 0), and x(v + c) - x(v) = dt/2 c
            space->create_vector(iterate);
            iterate.set(0.);
            fill_x_component(correction, 2 * c.gap / c.dt);
        } else {
            iterate = x0;
            fill_x_component(correction, c.gap);
        }

        Scalar_t alpha = -1;
        utopia_test_assert(fun->get_alpha(*fun, g, iterate, correction, alpha));
        return alpha;
    }

    // The predictor x_pred = x_n + dt v_hat_n, v_hat_n = T^-1 min(T v_n, gap / dt), of the stabilized schemes,
    // recovered from time_derivative(0) = v_hat_n - 2/dt x_pred = -(x_pred + x_n) / dt
    template <class Scheme>
    void predictor(Scheme &fun, const Scalar_t dt, Vector_t &x_pred) {
        Vector_t zero(layout(fun.velocity()), 0.), v;
        fun.time_derivative(zero, v);
        x_pred = -dt * v;
        x_pred -= fun.solution();
    }

    // After the step 0 -> x1 = 0.04 e_x, the velocity v1 = 2/dt x1 moves the right face by dt v1 = 0.08,
    // while only 0.01 is left to the obstacle (gap 0.05). The predictor must stop at the obstacle.
    template <class Scheme>
    void check_stabilized_predictor(const std::string &name, const bool velocity_unknown, const bool trivial) {
        Config c;
        c.gap = 0.05;
        c.with_bc = false;
        c.trivial_obstacle = trivial;

        auto space = make_space(c.with_bc);
        auto fun = make_scheme<Scheme>(space, c);

        Vector_t x0, x1, v1;
        space->create_vector(x0);
        x0.set(0.);
        utopia_test_assert(fun->setup_IVP(x0));

        space->create_vector(x1);
        fill_x_component(x1, 0.04);
        v1 = (2 / c.dt) * x1;

        utopia_test_assert(fun->update_IVP(velocity_unknown ? v1 : x1));

        Vector_t x_pred;
        predictor(*fun, c.dt, x_pred);

        // Expected predictor, with the gap measured in the configuration x1
        auto obstacle = make_obstacle(c.gap);
        space->displace(x1);
        obstacle->assemble(*space);
        space->displace(-x1);

        Vector_t t_inc, expected_inc;
        obstacle->transform(Vector_t(c.dt * v1), t_inc);
        t_inc = utopia::min(t_inc, obstacle->gap());
        obstacle->inverse_transform(t_inc, expected_inc);

        const Scalar_t v_error = norm_infty(Vector_t(fun->velocity() - v1));
        const Scalar_t pred_error = norm_infty(Vector_t(x_pred - x1 - expected_inc));
        const Scalar_t clamped = norm_infty(Vector_t(c.dt * v1 - expected_inc));
        const Scalar_t violation = constraint_violation(*space, *obstacle, x1, Vector_t(x_pred - x1));

        if (!space->comm().rank()) {
            utopia::out() << name << ": |v_1 - 2/dt x_1| = " << v_error << ", |x_pred - expected| = " << pred_error
                          << ", |dt v_1 - (x_pred - x_1)| = " << clamped << ", constraint violation = " << violation
                          << "\n";
        }

        utopia_test_assert(v_error <= 1e-12);
        utopia_test_assert(clamped >= 0.05);
        utopia_test_assert(pred_error <= 1e-12);
        utopia_test_assert(violation <= 1e-12);
    }

    // A body force pushes the cube into the obstacle (free flight would move it by 0.125 >> gap at t = 0.5).
    // Every step must be admissible, and the body must reach the obstacle. Same time loop as NLSolve.
    // The default settings of the schemes are used (in particular the line-search damping)
    template <class Scheme>
    void contact_run(const std::string &name) {
        Config c;
        c.gap = 0.02;
        c.n_steps = 16;
        c.body_force = 1;
        c.with_bc = false;
        c.enable_line_search = true;
        c.damping = 0;
        c.min_barrier_parameter = 1e-10;

        auto space = make_space(c.with_bc);
        auto fun = make_scheme<Scheme>(space, c);

        Vector_t x;
        initial_displacement(*space, x);
        utopia_test_assert(fun->setup_IVP(x));

        // The line search exists after setup_IVP (as in NLSolve)
        auto newton = make_newton();
        auto line_search = fun->line_search();
        utopia_test_assert(static_cast<bool>(line_search));
        newton->set_line_search_strategy(line_search);

        auto obstacle = make_obstacle(c.gap);
        Vector_t zero(layout(x), 0.);

        // Smallest distance gap - T x of the right face from the obstacle, over all time steps
        Scalar_t min_distance = c.gap;
        int n_steps = 0;

        do {
            fun->begin_time_step();
            space->apply_constraints_update(x);
            fun->initial_guess_for_solver(x);
            utopia_test_assert(newton->solve(*fun, x));
            utopia_test_assert(fun->update_IVP(x));

            const Scalar_t distance = -constraint_violation(*space, *obstacle, zero, fun->solution());
            min_distance = std::min(min_distance, distance);
            ++n_steps;
        } while (!fun->is_IVP_solved());

        if (!space->comm().rank()) {
            utopia::out() << name << ": steps = " << n_steps
                          << ", smallest distance to the obstacle = " << min_distance << " (gap " << c.gap << ")\n";
        }

        utopia_test_assert(n_steps == c.n_steps);
        utopia_test_assert(min_distance >= 0.);
        utopia_test_assert(min_distance < 0.25 * c.gap);
    }

    // Largest normal velocity |(T v)_i| over the contact dofs that are on the obstacle at x_old and at x_new, i.e.,
    // that did not move in normal direction. n_resting is the number of such dofs
    Scalar_t resting_normal_velocity(FS_t &space,
                                     Contact_t &obstacle,
                                     const Vector_t &x_old,
                                     const Vector_t &x_new,
                                     const Vector_t &v,
                                     Scalar_t &n_resting) {
        // Gap in the reference configuration
        Vector_t zero(layout(v), 0.);
        space.displace(zero);
        obstacle.assemble(space);

        Vector_t t_old, t_new, t_v;
        obstacle.transform(x_old, t_old);
        obstacle.transform(x_new, t_new);
        obstacle.transform(v, t_v);

        // 1 for the resting dofs and 0 elsewhere
        Vector_t resting(layout(v), 0.);
        {
            auto gap_view = const_local_view_device(obstacle.gap());
            auto contact_view = const_local_view_device(obstacle.is_contact());
            auto old_view = const_local_view_device(t_old);
            auto new_view = const_local_view_device(t_new);
            auto resting_view = local_view_device(resting);

            parallel_for(
                local_range_device(resting), UTOPIA_LAMBDA(const SizeType_t i) {
                    const Scalar_t g = gap_view.get(i);
                    const bool on_obstacle = contact_view.get(i) > 0.5 && g - old_view.get(i) <= 1e-12 &&
                                             g - new_view.get(i) <= 1e-12;
                    resting_view.set(i, on_obstacle ? 1. : 0.);
                });
        }

        n_resting = sum(resting);
        Vector_t resting_v = e_mul(resting, t_v);
        return norm_infty(resting_v);
    }

    // Total energy 1/2 v^T M v + 1/2 x^T K x - f^T x (linear elasticity and constant body force)
    template <class Scheme>
    Scalar_t total_energy(const Scheme &fun) {
        const Vector_t &x = fun.solution();
        const Vector_t &v = fun.velocity();

        Vector_t zero(layout(x), 0.), minus_f;
        Matrix_t K;
        utopia_test_assert(fun.function()->hessian(x, K));
        utopia_test_assert(fun.function()->gradient(zero, minus_f));

        Vector_t Mv = (*fun.mass_matrix()) * v;
        Vector_t Kx = K * x;
        const Scalar_t kinetic = 0.5 * dot(v, Mv);
        const Scalar_t elastic = 0.5 * dot(x, Kx);
        const Scalar_t load = dot(minus_f, x);
        return kinetic + elastic + load;
    }

    // Largest relative change of the total energy over the time steps. Without contact and with a constant
    // load the Newmark method (gamma = 1/2, beta = 1/4) conserves the energy exactly
    template <class Scheme>
    Scalar_t energy_drift(const Config &c) {
        auto space = make_space(c.with_bc);
        auto fun = make_scheme<Scheme>(space, c);
        auto newton = make_newton();

        Vector_t x;
        initial_displacement(*space, x);
        utopia_test_assert(fun->setup_IVP(x));

        const Scalar_t energy_0 = total_energy(*fun);
        Scalar_t drift = 0;

        do {
            fun->begin_time_step();
            space->apply_constraints_update(x);
            fun->initial_guess_for_solver(x);
            utopia_test_assert(newton->solve(*fun, x));
            utopia_test_assert(fun->update_IVP(x));
            drift = std::max(drift, std::abs(total_energy(*fun) - energy_0));
        } while (!fun->is_IVP_solved());

        return drift / std::abs(energy_0);
    }

}  // namespace

// Without contact (far obstacle, the bounded barrier vanishes) all six Newmark schemes solve the same
// discrete problem: same displacement, and the time-dependent Dirichlet data imposed at t_{n+1}
void moonolith_stk_newmark_schemes_agree() {
    Config c;

    Scalar_t bc_ref = 0, bc = 0;
    Vector_t x_ref = run<Newmark_t>(c, bc_ref);
    const Scalar_t x_norm = norm_infty(x_ref);
    utopia_test_assert(bc_ref <= 1e-12);

    auto check = [&](const std::string &name, const Vector_t &x, const Scalar_t bc_violation) {
        const Scalar_t diff = norm_infty(Vector_t(x - x_ref)) / x_norm;

        if (!x.comm().rank()) {
            utopia::out() << "moonolith_stk_newmark_schemes_agree: " << name
                          << ": |x - x_newmark| / |x_newmark| = " << diff << ", BC violation = " << bc_violation
                          << "\n";
        }

        utopia_test_assert(diff <= 1e-8);
        utopia_test_assert(bc_violation <= 1e-12);
    };

    check("VelocityNewmarkIntegrator", run<VelocityNewmark_t>(c, bc), bc);
    check("ObstacleNewmark", run<ObstacleNewmark_t>(c, bc), bc);
    check("ObstacleVelocityNewmark", run<ObstacleVelocityNewmark_t>(c, bc), bc);
    check("ObstacleStabilizedNewmark", run<StabilizedNewmark_t>(c, bc), bc);
    check("ObstacleStabilizedVelocityNewmark", run<StabilizedVelocityNewmark_t>(c, bc), bc);
}

// Without contact (far obstacle) every scheme conserves the energy exactly. This requires, among others,
// a symmetric mass matrix and consistent initial data
void moonolith_stk_newmark_energy() {
    Config c;
    c.n_steps = 16;
    c.body_force = 1;
    c.with_bc = false;

    auto check = [&](const std::string &name, const Scalar_t drift) {
        if (!mpi_world_rank()) {
            utopia::out() << "moonolith_stk_newmark_energy: " << name << ": max |E_n - E_0| / |E_0| = " << drift
                          << "\n";
        }

        utopia_test_assert(drift <= 1e-10);
    };

    check("NewmarkIntegrator", energy_drift<Newmark_t>(c));
    check("VelocityNewmarkIntegrator", energy_drift<VelocityNewmark_t>(c));
    check("ObstacleNewmark", energy_drift<ObstacleNewmark_t>(c));
    check("ObstacleVelocityNewmark", energy_drift<ObstacleVelocityNewmark_t>(c));
    check("ObstacleStabilizedNewmark", energy_drift<StabilizedNewmark_t>(c));
    check("ObstacleStabilizedVelocityNewmark", energy_drift<StabilizedVelocityNewmark_t>(c));
}

void moonolith_stk_newmark_obstacle_derivatives() {
    check_obstacle_derivatives<ObstacleNewmark_t>("ObstacleNewmark", 0.0005, 0.0005);
    check_obstacle_derivatives<ObstacleVelocityNewmark_t>("ObstacleVelocityNewmark", 0.01, 0.01);
    check_obstacle_derivatives<StabilizedVelocityNewmark_t>("ObstacleStabilizedVelocityNewmark", 0.01, 0.01);

    for (const std::string function_type : {"BoundedLogBarrier", "ShiftedPenalty"}) {
        check_obstacle_derivatives<ObstacleNewmark_t>("ObstacleNewmark", 0.0005, 0.0005, function_type, 0.1);
        check_obstacle_derivatives<ObstacleVelocityNewmark_t>("ObstacleVelocityNewmark", 0.01, 0.01, function_type, 0.1);
        check_obstacle_derivatives<StabilizedVelocityNewmark_t>(
            "ObstacleStabilizedVelocityNewmark", 0.01, 0.01, function_type, 0.1);
    }

    check_obstacle_derivatives<ObstacleVelocityNewmark_t>("ObstacleVelocityNewmark", 0.01, 0.01, "ShiftedPenalty");

    // Curved obstacle
    for (const std::string function_type : {"BoundedLogBarrier", "ShiftedPenalty"}) {
        check_obstacle_derivatives<ObstacleNewmark_t>("ObstacleNewmark", 0.0005, 0.0005, function_type, 0, 1);
        check_obstacle_derivatives<ObstacleVelocityNewmark_t>(
            "ObstacleVelocityNewmark", 0.01, 0.01, function_type, 0, 1);
        check_obstacle_derivatives<StabilizedVelocityNewmark_t>(
            "ObstacleStabilizedVelocityNewmark", 0.01, 0.01, function_type, 0, 1);
    }
}

void moonolith_stk_newmark_obstacle_line_search() {
    Config c;
    c.gap = 0.05;
    c.with_bc = false;
    c.enable_line_search = true;

    auto check = [&](const std::string &name, const Scalar_t alpha) {
        if (!mpi_world_rank()) {
            utopia::out() << "moonolith_stk_newmark_obstacle_line_search: " << name << ": alpha = " << alpha
                          << " (admissible 0.7)\n";
        }

        utopia_test_assert(alpha <= 0.7 + 1e-10);
        utopia_test_assert(alpha >= 0.6);
    };

    for (bool trivial : {false, true}) {
        c.trivial_obstacle = trivial;
        const std::string suffix = trivial ? " (trivial_obstacle)" : "";

        check("ObstacleNewmark" + suffix, line_search_alpha<ObstacleNewmark_t>(c, false));
        check("ObstacleVelocityNewmark" + suffix, line_search_alpha<ObstacleVelocityNewmark_t>(c, true));
        check("ObstacleStabilizedVelocityNewmark" + suffix, line_search_alpha<StabilizedVelocityNewmark_t>(c, true));
    }
}

void moonolith_stk_newmark_stabilized_predictor() {
    check_stabilized_predictor<StabilizedNewmark_t>("ObstacleStabilizedNewmark", false, false);
    check_stabilized_predictor<StabilizedVelocityNewmark_t>("ObstacleStabilizedVelocityNewmark", true, false);
    check_stabilized_predictor<StabilizedVelocityNewmark_t>(
        "ObstacleStabilizedVelocityNewmark (trivial_obstacle)", true, true);
}

void moonolith_stk_newmark_contact() {
    contact_run<ObstacleNewmark_t>("moonolith_stk_newmark_contact: ObstacleNewmark");
    contact_run<ObstacleVelocityNewmark_t>("moonolith_stk_newmark_contact: ObstacleVelocityNewmark");
    contact_run<StabilizedVelocityNewmark_t>("moonolith_stk_newmark_contact: ObstacleStabilizedVelocityNewmark");
}

UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_schemes_agree);
UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_energy);
UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_obstacle_derivatives);
UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_obstacle_line_search);
UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_stabilized_predictor);
// ObstacleStabilizedNewmark with the box-constrained SQP solver and ObstacleFEFunction (as in ObstacleApp).
// A body force pushes the cube into the obstacle: every step must be admissible, the body must reach the
// obstacle, and the energy must not increase (it is conserved without contact and dissipated at contact).
// Nodes that rest on the obstacle must have zero normal velocity (with the velocity update of Algorithm (21) in
// Deuflhard, Krause, Ertel, IJNME 2008, they keep the impact velocity)
void moonolith_stk_newmark_stabilized_contact() {
    Config c;
    c.gap = 0.02;
    c.n_steps = 16;
    c.body_force = 1;
    c.with_bc = false;

    auto space = make_space(c.with_bc);
    auto fun = std::make_shared<StabilizedNewmark_t>(std::make_shared<FEModelFunction<FS_t>>(space));
    auto obs_fun = std::make_shared<ObstacleFEFunction<FS_t>>(fun);

    // Reads the time integrator and the obstacle, which is passed to the time integrator
    auto params = make_params(c);
    obs_fun->read(params);

    BoxConstrainedFEFunctionSolver<FS_t> solver;
    solver.verbose(false);
    auto qp = std::make_shared<MPRGP<Matrix_t, Vector_t>>();
    qp->atol(1e-14);
    qp->rtol(1e-14);
    qp->stol(1e-14);
    qp->max_it(10000);
    solver.qp_solver_ = qp;

    Vector_t x;
    initial_displacement(*space, x);
    utopia_test_assert(fun->setup_IVP(x));

    auto obstacle = make_obstacle(c.gap);
    Vector_t zero(layout(x), 0.);

    const Scalar_t energy_0 = total_energy(*fun);
    Scalar_t energy_prev = energy_0;
    Scalar_t max_energy_increase = 0;
    Scalar_t min_distance = c.gap;
    Scalar_t max_resting_velocity = 0, n_resting = 0;
    int n_steps = 0;

    do {
        fun->begin_time_step();
        const Vector_t x_old = fun->solution();
        utopia_test_assert(solver.solve(*obs_fun, x));
        utopia_test_assert(fun->update_IVP(x));

        Scalar_t n_resting_step = 0;
        max_resting_velocity = std::max(
            max_resting_velocity,
            resting_normal_velocity(*space, *obstacle, x_old, fun->solution(), fun->velocity(), n_resting_step));
        n_resting += n_resting_step;

        const Scalar_t distance = -constraint_violation(*space, *obstacle, zero, fun->solution());
        min_distance = std::min(min_distance, distance);

        const Scalar_t energy = total_energy(*fun);
        max_energy_increase = std::max(max_energy_increase, energy - energy_prev);
        energy_prev = energy;
        ++n_steps;
    } while (!fun->is_IVP_solved());

    if (!space->comm().rank()) {
        utopia::out() << "moonolith_stk_newmark_stabilized_contact: steps = " << n_steps
                      << ", smallest distance to the obstacle = " << min_distance << " (gap " << c.gap
                      << "), energy: initial = " << energy_0 << ", final = " << energy_prev
                      << ", largest increase in a step = " << max_energy_increase
                      << ", largest normal velocity of resting nodes = " << max_resting_velocity << " ("
                      << n_resting << " resting dofs over all steps)\n";
    }

    utopia_test_assert(n_resting > 0);
    utopia_test_assert(max_resting_velocity <= 1e-10);

    utopia_test_assert(n_steps == c.n_steps);
    utopia_test_assert(min_distance >= -1e-10);
    utopia_test_assert(min_distance < 1e-6);
    utopia_test_assert(max_energy_increase <= 1e-8 * std::abs(energy_0));
    utopia_test_assert(energy_prev < energy_0);
}

UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_contact);
UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_stabilized_contact);

namespace {

    // Runs the scheme with NLSolve, which writes the output at every step, and compares the velocity field
    // of the last step on disk with the velocity of the scheme
    template <class Scheme>
    void check_velocity_output(const std::string &name) {
        Config c;
        c.n_steps = 4;
        const std::string path = "newmark_velocity_output_test.e";

        auto space = make_space(c.with_bc);
        Vector_t velocity;

        {
            auto fun = std::make_shared<Scheme>(std::make_shared<FEModelFunction<FS_t>>(space));
            auto params = make_params(c);
            params.set("output_path", path);
            fun->read(params);
            attach_obstacle(*fun, c);

            auto x0 = std::make_shared<Vector_t>();
            initial_displacement(*space, *x0);

            NLSolve<FS_t> nl;
            nl.verbose(false);
            nl.init(fun);
            nl.set_solution(std::make_shared<Field<FS_t>>("x", space, x0));
            nl.set_solver(make_newton());
            utopia_test_assert(nl.solve());

            velocity = fun->velocity();
        }

        IO<FS_t> input(*space);
        input.set_input_path(path);
        input.import_all_field_data(true);
        utopia_test_assert(input.open_input());
        utopia_test_assert(input.load_time_step(c.dt * c.n_steps));

        Field<FS_t> velocity_disk("velocity", space);
        utopia_test_assert(input.read_nodal(velocity_disk));

        const Scalar_t diff = norm_infty(Vector_t(velocity_disk.data() - velocity));
        const Scalar_t v_norm = norm_infty(velocity);

        if (!space->comm().rank()) {
            utopia::out() << "moonolith_stk_newmark_stabilized_output: " << name << ", |v| = " << v_norm
                          << ", |v_disk - v| = " << diff << "\n";
        }

        utopia_test_assert(v_norm > 0.);
        utopia_test_assert(diff <= 1e-12 * v_norm);
    }

}  // namespace

// The stabilized schemes write the velocity field
void moonolith_stk_newmark_stabilized_output() {
    check_velocity_output<StabilizedNewmark_t>("ObstacleStabilizedNewmark");
    check_velocity_output<StabilizedVelocityNewmark_t>("ObstacleStabilizedVelocityNewmark");
}

UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_stabilized_output);

namespace {

    struct Trajectory {
        Vector_t x, v;
        std::vector<Scalar_t> energy;
    };

    Config csn_config() {
        Config c;
        c.gap = 0.02;
        c.n_steps = 16;
        c.body_force = 1;
        c.with_bc = false;
        return c;
    }

    // Contact-stabilized Newmark with the exact constraint (ObstacleStabilizedNewmark, ObstacleFEFunction, MPRGP)
    Trajectory csn_exact(const Config &c) {
        auto space = make_space(c.with_bc);
        auto fun = std::make_shared<StabilizedNewmark_t>(std::make_shared<FEModelFunction<FS_t>>(space));
        auto obs_fun = std::make_shared<ObstacleFEFunction<FS_t>>(fun);
        auto params = make_params(c);
        obs_fun->read(params);

        BoxConstrainedFEFunctionSolver<FS_t> solver;
        solver.verbose(false);
        auto qp = std::make_shared<MPRGP<Matrix_t, Vector_t>>();
        qp->atol(1e-14);
        qp->rtol(1e-14);
        qp->stol(1e-14);
        qp->max_it(10000);
        solver.qp_solver_ = qp;

        Vector_t x;
        initial_displacement(*space, x);
        utopia_test_assert(fun->setup_IVP(x));

        Trajectory t;
        t.energy.push_back(total_energy(*fun));
        do {
            fun->begin_time_step();
            utopia_test_assert(solver.solve(*obs_fun, x));
            utopia_test_assert(fun->update_IVP(x));
            t.energy.push_back(total_energy(*fun));
        } while (!fun->is_IVP_solved());

        t.x = fun->solution();
        t.v = fun->velocity();
        return t;
    }

    // Same problem with ObstacleStabilizedVelocityNewmark (barrier with parameter mu)
    Trajectory csn_barrier(Config c, const Scalar_t mu, const Scalar_t min_mu) {
        c.barrier_parameter = mu;
        c.enable_line_search = true;

        auto space = make_space(c.with_bc);
        auto fun = std::make_shared<StabilizedVelocityNewmark_t>(std::make_shared<FEModelFunction<FS_t>>(space));
        auto params = make_params(c);
        params.set("min_barrier_parameter", min_mu);
        fun->read(params);
        Vector_t x;
        initial_displacement(*space, x);
        utopia_test_assert(fun->setup_IVP(x));

        // The line search exists after setup_IVP (as in NLSolve)
        auto newton = make_newton();
        auto line_search = fun->line_search();
        utopia_test_assert(static_cast<bool>(line_search));
        newton->set_line_search_strategy(line_search);

        Trajectory t;
        t.energy.push_back(total_energy(*fun));
        do {
            fun->begin_time_step();
            fun->initial_guess_for_solver(x);
            utopia_test_assert(newton->solve(*fun, x));
            utopia_test_assert(fun->update_IVP(x));
            t.energy.push_back(total_energy(*fun));
        } while (!fun->is_IVP_solved());

        t.x = fun->solution();
        t.v = fun->velocity();
        return t;
    }

}  // namespace

// ObstacleStabilizedVelocityNewmark is the contact-stabilized Newmark method (Deuflhard, Krause, Ertel) with the
// contact constraint replaced by a log-barrier. Its trajectory must converge to the one of the exact method
// (ObstacleStabilizedNewmark with the box-constrained solver) as the barrier parameter mu goes to zero
void moonolith_stk_newmark_stabilized_velocity_vs_csn() {
    const Config c = csn_config();
    const Trajectory ref = csn_exact(c);

    auto distance = [&](const Trajectory &t, Scalar_t &ex, Scalar_t &ev) {
        ex = Scalar_t(norm_infty(Vector_t(t.x - ref.x))) / Scalar_t(norm_infty(ref.x));
        ev = Scalar_t(norm_infty(Vector_t(t.v - ref.v))) / Scalar_t(norm_infty(ref.v));
    };

    // Fixed barrier parameter: the displacement error is O(mu)
    Scalar_t ex_prev = 0;
    for (Scalar_t mu : {1e-5, 1e-6, 1e-7}) {
        Scalar_t ex, ev;
        distance(csn_barrier(c, mu, mu), ex, ev);

        if (!mpi_world_rank()) {
            utopia::out() << "moonolith_stk_newmark_stabilized_velocity_vs_csn: mu = " << mu
                          << ", |x - x_csn| / |x_csn| = " << ex << ", |v - v_csn| / |v_csn| = " << ev << "\n";
        }

        if (ex_prev > 0) {
            utopia_test_assert(ex < 0.2 * ex_prev);
        }

        utopia_test_assert(ev < 1e-3);

        ex_prev = ex;
    }

    utopia_test_assert(ex_prev < 1e-4);

    // Default usage: the barrier parameter shrinks during the Newton iterations (1e-3 -> 1e-10)
    Scalar_t ex, ev;
    distance(csn_barrier(c, 1e-3, 1e-10), ex, ev);

    if (!mpi_world_rank()) {
        utopia::out() << "moonolith_stk_newmark_stabilized_velocity_vs_csn: mu = 1e-3 -> 1e-10"
                      << ", |x - x_csn| / |x_csn| = " << ex << ", |v - v_csn| / |v_csn| = " << ev << "\n";
    }

    utopia_test_assert(ex < 1e-4);
    utopia_test_assert(ev < 1e-4);
}

UTOPIA_REGISTER_TEST_FUNCTION(moonolith_stk_newmark_stabilized_velocity_vs_csn);

#endif  // UTOPIA_ENABLE_INTREPID2 && UTOPIA_ENABLE_TINY_EXPR
