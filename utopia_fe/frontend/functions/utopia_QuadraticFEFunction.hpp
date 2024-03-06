#ifndef UTOPIA_QUADRATIC_FE_FUNCTION_HPP
#define UTOPIA_QUADRATIC_FE_FUNCTION_HPP

namespace utopia {

    template <class FunctionSpace>
    class QuadraticFEFunction : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;

        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using SimulationTime = utopia::SimulationTime<Scalar_t>;

        QuadraticFEFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &function) : function_(function) {
            assert(function_);
            H0 = utopia::make_unique<Matrix_t>();
            g0 = utopia::make_unique<Vector_t>();
        }

        inline Communicator_t &comm() override { return H0->comm(); }

        inline const Communicator_t &comm() const override { return H0->comm(); }

        void read(Input &in) override {
            function_->read(in);
            in.get("export_tensors", export_tensors_);
        }

        bool value(const Vector_t &x, Scalar_t &v) const override {
            Vector_t g = (*H0) * x;
            v = dot(x, g) + dot(x, (*g0));
            return true;
        }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            if (!initialize()) return false;

            g = (*H0) * x;
            g += (*g0);

            if (must_apply_constraints_) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        bool update(const Vector_t &x) override { return function_->update(x); }

        bool hessian(const Vector_t &, Matrix_t &H) const override {
            if (!initialize()) return false;

            H = (*H0);
            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            if (!initialize()) return false;

            H = (*H0);
            g = H * x;
            g += (*g0);

            if (must_apply_constraints_) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        void create_solution_vector(Vector_t &x) override { function_->create_solution_vector(x); }

        void apply_constraints(Vector_t &x) const override { function_->apply_constraints(x); }

        void set_environment(const std::shared_ptr<Environment_t> &env) override { function_->set_environment(env); }

        const std::shared_ptr<Matrix_t> &mass_matrix() const override { return function_->mass_matrix(); }

        bool assemble_mass_matrix() override { return function_->assemble_mass_matrix(); }

        bool assemble_mass_matrix(Matrix_t &mass_matrix) override {
            return function_->assemble_mass_matrix(mass_matrix);
        }

        const std::shared_ptr<FunctionSpace> &space() const override { return function_->space(); }

        bool initialize_hessian(Matrix_t &H, Matrix_t &H_preconditioner) const override {
            return function_->initialize_hessian(H, H_preconditioner);
        }

        bool is_time_dependent() const override { return function_->is_time_dependent(); }

        bool is_linear() const override {
            assert(function_->is_linear());
            return true;
        }

        void must_apply_constraints_to_assembled(const bool value) override {
            function_->must_apply_constraints_to_assembled(value);
            must_apply_constraints_ = value;
        }

        bool report_solution(const Vector_t &x) override { return function_->report_solution(x); }

        bool update_IVP(const Vector_t &x) override { return function_->update_IVP(x); }
        bool setup_IVP(Vector_t &x) override { return function_->setup_IVP(x); }
        bool is_IVP_solved() override { return function_->is_IVP_solved(); }

        void set_time(const std::shared_ptr<SimulationTime> &time) override {
            assert(function_);
            if (function_) {
                function_->set_time(time);
            }
        }

        void post_solve(Vector_t &x) override {
            if (function_) function_->post_solve(x);
        }

    public:
        std::shared_ptr<FEFunctionInterface<FunctionSpace>> function_;

        std::unique_ptr<Matrix_t> H0;
        std::unique_ptr<Vector_t> g0;

        bool must_apply_constraints_{true};
        bool export_tensors_{false};

        bool initialize() const {
            if (!H0->empty()) {
                return true;
            }

            Vector_t x;
            function_->create_solution_vector(x);
            function_->create_solution_vector(*g0);

            x.set(0.0);
            g0->set(0.0);

            if (!function_->hessian_and_gradient(x, *H0, *g0)) {
                assert(false);
                return false;
            }

            rename("H0", *H0);
            rename("g0", *g0);

            if (export_tensors_) {
                write("load_H0.m", *H0);
                write("load_g0.m", *g0);
            }

            return true;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_QUADRATIC_FE_FUNCTION_HPP