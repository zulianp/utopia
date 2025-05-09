#ifndef UTOPIA_BOX_CONSTRAINED_FE_FUNCTION_HPP
#define UTOPIA_BOX_CONSTRAINED_FE_FUNCTION_HPP

#include "utopia_BoxConstraints.hpp"

#include "utopia_FEModelFunction.hpp"

#include "utopia_polymorphic_QPSolver.hpp"

#include "utopia_MonotoneSemiGeometricMultigrid.hpp"

#include "utopia_polymorphic_QPSolver_impl.hpp"

namespace utopia {

    template <class FunctionSpace>
    class BoxConstrainedFEFunction : public FEFunctionInterface<FunctionSpace> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using SimulationTime = utopia::SimulationTime<Scalar_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;

        ~BoxConstrainedFEFunction() override = default;

        virtual bool has_nonlinear_constraints() const = 0;
        virtual bool constraints_gradient(const Vector_t &x, BoxConstraints<Vector_t> &box) = 0;
        virtual bool has_transformation() const = 0;
        virtual void transform(const Vector_t &x, Vector_t &x_constrained) = 0;
        virtual void inverse_transform(const Vector_t &x_constrained, Vector_t &x) = 0;
        virtual void transform(const Matrix_t &H, Matrix_t &H_constrained) = 0;
        virtual std::shared_ptr<Vector_t> selection() = 0;

        bool value(const Vector_t & /*point*/, Scalar_t & /*value*/) const override { return false; }

        bool gradient(const Vector_t &x, Vector_t &g) const override { return unconstrained_->gradient(x, g); }
        bool hessian(const Vector_t &x, Matrix_t &H) const override { return unconstrained_->hessian(x, H); }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            return unconstrained_->hessian_and_gradient(x, H, g);
        }

        void read(Input &in) override { unconstrained_->read(in); }

        void create_solution_vector(Vector_t &x) override { return unconstrained_->create_solution_vector(x); }

        void apply_constraints(Vector_t &x) const override { return unconstrained_->apply_constraints(x); }
        void apply_constraints(const Vector_t &u, const Scalar_t scale_factor, Vector_t &in_out) {
            return unconstrained()->space().apply_constraints(u, scale_factor, in_out);
        }

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            unconstrained_->set_environment(env);
        }

        const std::shared_ptr<Matrix_t> &mass_matrix() const override { return unconstrained_->mass_matrix(); }
        bool assemble_mass_matrix() override { return unconstrained_->assemble_mass_matrix(); }

        bool assemble_mass_matrix(Matrix_t &mass_matrix) override {
            return unconstrained_->assemble_mass_matrix(mass_matrix);
        }

        const std::shared_ptr<FunctionSpace> &space() const override { return unconstrained_->space(); }

        bool is_time_dependent() const override { return unconstrained_->is_time_dependent(); }
        bool is_linear() const override { return unconstrained_->is_linear(); }

        void must_apply_constraints_to_assembled(const bool val) override {
            unconstrained_->must_apply_constraints_to_assembled(val);
        }

        bool report_solution(const Vector_t &x) override { return unconstrained_->report_solution(x); }

        bool update_IVP(const Vector_t &x) override { return unconstrained_->update_IVP(x); }
        bool setup_IVP(Vector_t &x) override { return unconstrained_->setup_IVP(x); }
        bool setup_IVP(IO<FunctionSpace> &input) override { return unconstrained_->setup_IVP(input); }
        bool is_IVP_solved() override { return unconstrained_->is_IVP_solved(); }

        void set_output_db(const std::shared_ptr<IO<FunctionSpace>> &io) override {
            assert(unconstrained_);

            if (!unconstrained_) {
                Utopia::Abort("BoxConstrainedFEFunction::set_output_db called on uninitialized object!");
            }

            unconstrained_->set_output_db(io);
        }

        bool register_output(IO<FunctionSpace> &io) override {
            assert(unconstrained_);

            if (!unconstrained_) {
                Utopia::Abort("BoxConstrainedFEFunction::register_output called on uninitialized object!");
            }

            return unconstrained_->register_output(io);
        }

        bool update_output(IO<FunctionSpace> &io) override {
            assert(unconstrained_);

            if (!unconstrained_) {
                Utopia::Abort("BoxConstrainedFEFunction::update_output called on uninitialized object!");
            }

            return unconstrained_->update_output(io);
        }

        BoxConstrainedFEFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained) {
            initialize(unconstrained);
        }

        BoxConstrainedFEFunction() {}

        virtual void initialize(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained) {
            unconstrained_ = unconstrained;
            assert(unconstrained_);
        }

        virtual bool is_initialized() const { return static_cast<bool>(unconstrained_); }

        void set_time(const std::shared_ptr<SimulationTime> &time) override {
            assert(unconstrained_);
            if (unconstrained_) {
                unconstrained_->set_time(time);
            }
        }

        inline Communicator_t &comm() override { return unconstrained_->comm(); }

        inline const Communicator_t &comm() const override { return unconstrained_->comm(); }

        const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained() const { return unconstrained_; };

    private:
        std::shared_ptr<FEFunctionInterface<FunctionSpace>> unconstrained_;
    };

    template <class FunctionSpace>
    class BoxConstrainedFEFunctionSolver final : public NonLinearSolver<typename Traits<FunctionSpace>::Vector> {
    public:
        using Super = utopia::NonLinearSolver<typename Traits<FunctionSpace>::Vector>;
        using Function_t = utopia::BoxConstrainedFEFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using QPSolver_t = utopia::QPSolver<Matrix_t, Vector_t>;
        using OmniQPSolver_t = utopia::OmniQPSolver<Matrix_t, Vector_t>;

        void read(Input &in) override {
            Super::read(in);

            ensure_qp_solver();

            in.get("qp_solver", *qp_solver_);
            in.get("material_iter_tol", material_iter_tol_);
            in.get("max_constraints_iterations", max_constraints_iterations_);
            in.get("rescale", rescale_);
            in.get("inverse_diagonal_scaling", inverse_diagonal_scaling_);
            in.get("print_active_set", print_active_set_);
            in.get("damping", damping_);
            in.get("export_sqp_system", export_sqp_system_);
        }

        void ensure_qp_solver() {
            if (!qp_solver_) {
                auto omni = std::make_shared<OmniQPSolver_t>();
                qp_solver_ = omni;
            }
        }

        bool solve(Function_t &fun, Vector_t &x) {
            UTOPIA_TRACE_SCOPE("BoxConstrainedFEFunctionSolver::solve");

            ensure_qp_solver();

            BoxConstraints<Vector_t> box;

            int max_material_iterations = this->max_it();

            if (fun.is_linear()) {
                max_material_iterations = 1;
            }

            Matrix_t H;
            Vector_t g, increment, x_old;

            // In case there is a system transformation to be performed
            Matrix_t H_c;
            Vector_t g_c, increment_c;

            fun.space()->create_matrix(H);
            fun.space()->create_vector(g);
            fun.space()->create_vector(increment);

            this->init_solver("BoxConstrainedFEFunctionSolver",
                              {" it. ", "|| u_old - u_new ||_2", "|| u_old - u_new ||_A"});

            bool first = true;
            int total_iter = 0;
            bool converged = false;
            for (int constraints_iter = 1; constraints_iter <= max_constraints_iterations_; ++constraints_iter) {
                x_old = x;
                fun.constraints_gradient(x, box);

                bool material_converged = false;
                converged = false;
                for (int material_iter = 1; material_iter <= max_material_iterations; ++material_iter, ++total_iter) {
                    fun.hessian_and_gradient(x, H, g);

                    if (g.has_nan_or_inf()) {
                        fun.space()->write("NaN.e", x);
                        Utopia::Abort("BoxConstrainedFEFunction: NaN found in gradient!\n");
                    }

                    if (rescale_ != 1.) {
                        // x.comm().root_print("Rescaling system with " + std::to_string(rescale_) + "\n");
                        H *= rescale_;
                        g *= rescale_;
                    }

                    // Use negative gradient instead
                    g *= -1;

                    fun.space()->apply_constraints(x, 1, g);
                    fun.space()->apply_constraints(H);

                    if (first) {
                        // fun.space()->apply_constraints(H, g);
                        first = false;
                    }
                    // else {
                    //     fun.space()->apply_constraints(H);
                    //     fun.space()->apply_zero_constraints(g);
                    // }

                    qp_solver_->set_box_constraints(box);

                    auto selection = fun.selection();
                    if (selection) {
                        qp_solver_->set_selection(selection);
                    }

                    bool qp_solver_converged = false;

                    if (fun.has_transformation()) {
                        fun.transform(H, H_c);
                        fun.transform(g, g_c);

                        if (empty(increment_c)) {
                            increment_c.zeros(layout(g_c));
                        } else {
                            increment_c.set(0.);
                        }

                        if (inverse_diagonal_scaling_) {
                            Vector_t d = diag(H_c);
                            d = 1. / d;
                            H_c.diag_scale_left(d);
                            g_c = e_mul(g_c, d);
                        }

                        if (export_sqp_system_) {
                            static int sys_num = 0;
                            rename("A", H_c);
                            rename("b", g_c);

                            Path dir("sqp");
                            if (!dir.exists()) {
                                dir.make_dir();
                            }

                            Path step_dir = dir / std::to_string(sys_num);

                            if (!step_dir.exists()) {
                                step_dir.make_dir();
                            }

                            write(step_dir / "load_A.m", H_c);
                            write(step_dir / "load_b.m", g_c);

                            if (box.upper_bound()) {
                                rename("ub", *box.upper_bound());
                                write(step_dir / "load_ub.m", *box.upper_bound());
                            }

                            if (box.lower_bound()) {
                                rename("lb", *box.lower_bound());
                                write(step_dir / "load_lb.m", *box.lower_bound());
                            }

                            sys_num++;
                        }

                        qp_solver_converged = qp_solver_->solve(H_c, g_c, increment_c);

                        if (damping_ != 1) {
                            increment_c *= damping_;
                        }

                        if (print_active_set_) {
                            // Count active nodes
                            auto count_a = box.count_active(increment_c, 1e-16);
                            if (increment.comm().rank() == 0) {
                                utopia::out() << "Active dofs: " << count_a << "\n";
                            }
                        }

                        fun.inverse_transform(increment_c, increment);

                        if (box.upper_bound()) {
                            *box.upper_bound() -= increment_c;
                        }

                        if (box.lower_bound()) {
                            *box.lower_bound() -= increment_c;
                        }

                    } else {
                        increment.set(0.0);

                        if (inverse_diagonal_scaling_) {
                            Vector_t d = diag(H);
                            d = 1. / d;
                            H.diag_scale_left(d);
                            g = e_mul(g, d);
                        }

                        qp_solver_converged = qp_solver_->solve(H, g, increment);

                        if (print_active_set_) {
                            // Count active nodes
                            auto count_a = box.count_active(increment, 1e-16);
                            if (increment.comm().rank() == 0) {
                                utopia::out() << "Active dofs: " << count_a << "\n";
                            }
                        }

                        if (box.upper_bound()) {
                            *box.upper_bound() -= increment;
                        }

                        if (box.lower_bound()) {
                            *box.lower_bound() -= increment;
                        }
                    }

                    x += increment;

                    if (!qp_solver_converged) {
                        assert(false);
                        utopia::err() << "BoxConstrainedFEFunctionSolver[Error] unable to solve QP problem "
                                         "terminating solution process\n";
                        return false;
                    }

                    if (fun.is_linear()) {
                        material_converged = qp_solver_converged;
                    } else {
                        Scalar_t material_inc_norm_A = std::sqrt(dot(increment, H * increment));

                        // Make sure that we do not consider the damping params when checking for convergence
                        material_inc_norm_A /= damping_;

                        if (this->verbose()) {
                            const Scalar_t material_inc_norm_2 = norm2(increment);

                            PrintInfo::print_iter_status(material_iter, {material_inc_norm_2, material_inc_norm_A});
                        }

                        if (material_inc_norm_A < material_iter_tol_) {
                            material_converged = true;
                            break;
                        }
                    }
                }

                Scalar_t x_diff_norm_A = std::sqrt(dot((x_old - x), H * (x_old - x)));

                if (this->verbose()) {
                    if (increment.comm().rank() == 0) {
                        utopia::out() << "Contact linearization (||u_old - u||_2):\n";
                    }

                    Scalar_t x_diff_norm_2 = norm2(x_old - x);
                    PrintInfo::print_iter_status(constraints_iter, {x_diff_norm_2, x_diff_norm_A});
                }

                converged = this->check_convergence(constraints_iter, 1, 1, x_diff_norm_A);
                if (converged) {
                    // Consider material convergence for overall convergence check
                    converged = material_converged;
                }

                if (converged && this->verbose()) {
                    x.comm().root_print("Converged!");
                    break;
                }

                if (!fun.has_nonlinear_constraints()) {
                    converged = material_converged;
                    break;
                }
            }

            return converged;
        }

        std::shared_ptr<QPSolver_t> qp_solver_;

        BoxConstrainedFEFunctionSolver() { register_fe_solvers(); }

    public:
        int max_constraints_iterations_{10};
        Scalar_t damping_{1};
        Scalar_t material_iter_tol_{1e-6};
        Scalar_t rescale_{1};
        bool inverse_diagonal_scaling_{false};
        bool print_active_set_{false};
        bool export_sqp_system_{false};

        // FIXME move somewhere else
        static void register_fe_solvers() {
            using MonotoneSemiGeometricMultigrid = utopia::MonotoneSemiGeometricMultigrid<FunctionSpace>;
            static bool registered = false;

            if (!registered) {
                registered = true;
                OmniQPSolver_t::registry().template register_solver<MonotoneSemiGeometricMultigrid>();
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_BOX_CONSTRAINED_FE_FUNCTION_HPP