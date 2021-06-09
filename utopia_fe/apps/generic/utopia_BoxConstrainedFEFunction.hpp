#ifndef UTOPIA_BOX_CONSTRAINED_FE_FUNCTION_HPP
#define UTOPIA_BOX_CONSTRAINED_FE_FUNCTION_HPP

#include "utopia_BoxConstraints.hpp"

#include "utopia_FEModelFunction.hpp"

#include "utopia_polymorphic_QPSolver.hpp"

namespace utopia {

    template <class FunctionSpace>
    class BoxConstrainedFEFunction : public FEFunctionInterface<FunctionSpace> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Environment_t = utopia::Environment<FunctionSpace>;

        ~BoxConstrainedFEFunction() override = default;

        virtual bool has_nonlinear_constraints() const = 0;
        virtual bool update_constraints(const Vector_t &x) = 0;
        virtual bool constraints_gradient(const Vector_t &x, BoxConstraints<Vector_t> &box) = 0;
        virtual bool has_transformation() const = 0;
        virtual void transform(const Vector_t &x, Vector_t &x_constrained) = 0;
        virtual void inverse_transform(const Vector_t &x_constrained, Vector_t &x) = 0;
        virtual void transform(const Matrix_t &H, Matrix_t &H_constrained) = 0;

        void read(Input &in) override { unconstrained_->read(in); }

        void create_solution_vector(Vector_t &x) override { return unconstrained_->create_solution_vector(x); }

        void apply_constraints(Vector_t &x) const override { return unconstrained_->apply_constraints(x); }

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
        bool is_IVP_solved() override { return unconstrained_->is_IVP_solved(); }

        BoxConstrainedFEFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
            : unconstrained_(unconstrained) {
            assert(unconstrained_);
        }

    private:
        std::shared_ptr<FEFunctionInterface<FunctionSpace>> unconstrained_;
    };

    template <class FunctionSpace>
    class BoxConstrainedFEFunctionSolver final : public NonLinearSolver<typename Traits<FunctionSpace>::Vector> {
    public:
        using Function_t = utopia::BoxConstrainedFEFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using QPSolver_t = utopia::QPSolver<Matrix_t, Vector_t>;
        using OmniQPSolver_t = utopia::OmniQPSolver<Matrix_t, Vector_t>;

        void read(Input &in) override {
            if (!qp_solver_) {
                auto omni = std::make_shared<OmniQPSolver_t>();
                qp_solver_ = omni;
            }

            in.get("qp_solver", *qp_solver_);
            in.get("update_factor", update_factor_);
            in.get("material_iter_tol", material_iter_tol_);
        }

        bool solve(Function_t &fun, Vector_t &x) {
            BoxConstraints<Vector_t> box;
            box.fill_empty_bounds(layout(x));
            fun.update_constraints(x);

            int max_material_iterations = this->max_it();

            if (fun.is_linear()) {
                max_material_iterations = 1;
            }

            Matrix_t H;
            Vector_t g, increment, x_old;

            // In case there is a system transformation to be performed
            Matrix_t H_c;
            Vector_t g_c, increment_c;

            this->space()->create_matrix(H);
            this->space()->create_vector(g);
            this->space()->create_vector(increment);

            this->init_solver("BoxConstrainedFEFunctionSolver", {" it. ", "|| u_old - u_new ||"});

            bool first = true;
            int total_iter = 0;
            for (int constraints_iter = 0; constraints_iter < max_constraints_iterations_; ++constraints_iter) {
                x_old = x;
                for (int material_iter = 0; material_iter < max_material_iterations; ++material_iter, ++total_iter) {
                    fun.hessian_and_gradient(x, H, g);
                    fun.constraints_gradient(x, box);

                    // Use negative gradient instead
                    g *= -1;

                    if (first) {
                        this->space()->apply_constraints(H, g);
                        first = false;
                    } else {
                        this->space()->apply_constraints(H);
                        this->space()->apply_zero_constraints(g);
                    }

                    qp_solver_->set_box_constraints(box);

                    if (fun.has_transformation()) {
                        fun.transform(H, H_c);
                        fun.transform(g, g_c);

                        if (empty(increment_c)) {
                            increment_c.zeros(layout(g_c));
                        } else {
                            increment_c.set(0.);
                        }

                        if (!qp_solver_->solve(H_c, g_c, increment_c)) {
                            assert(false);
                            utopia::err() << "BoxConstrainedFEFunctionSolver[Error] unable to solve QP problem "
                                             "terminating solution process\n";
                            return false;
                        }

                        fun.inverse_transform(increment_c, increment);

                    } else {
                        increment.set(0.0);
                        qp_solver_->solve(H, g, increment);
                    }

                    increment *= update_factor_;
                    x += increment;

                    const Scalar_t material_inc_norm = norm2(increment);

                    if (this->verbose()) {
                        PrintInfo::print_iter_status(total_iter, {material_inc_norm});
                    }

                    if (material_inc_norm < material_iter_tol_) {
                        break;
                    }
                }

                Scalar_t x_diff_norm = norm2(x_old - x);

                if (this->verbose()) {
                    PrintInfo::print_iter_status(total_iter, {x_diff_norm});
                }

                if (!fun.has_nonlinear_constraints()) {
                    break;
                }

                // Update constraints
                fun.update_constraints(x);
                fun.constraints_gradient(box);
            }
        }

        std::shared_ptr<QPSolver_t> qp_solver_;

    public:
        int max_constraints_iterations_{4};
        Scalar_t constraint_step_tol_{1e-4};
        Scalar_t update_factor_{1};
        Scalar_t material_iter_tol_{1e-6};
    };

}  // namespace utopia

#endif  // UTOPIA_BOX_CONSTRAINED_FE_FUNCTION_HPP