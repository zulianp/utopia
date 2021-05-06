#ifndef UTOPIA_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_NEWMARK_INTEGRATOR_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_MatrixTransformer.hpp"
#include "utopia_ProblemBase.hpp"

namespace utopia {

    template <class FunctionSpace>
    class TimeDependentFunction
        : public Function<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;

        virtual ~TimeDependentFunction() = default;

        virtual void set_environment(const std::shared_ptr<Environment_t> &env) = 0;

        virtual bool update_IVP(const Vector_t &x) = 0;
        virtual bool init_IVP(Vector_t &x) = 0;

        void read(Input &in) override {
            in.get("verbose", verbose_);
            in.get("delta_time", delta_time_);
        }

        inline Scalar_t delta_time() const { return delta_time_; }
        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) const { verbose_ = val; }

    private:
        Scalar_t delta_time_{0.1};
        bool verbose_{false};
    };

    template <class FunctionSpace>
    class NewmarkIntegrator final : public TimeDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using IO_t = utopia::IO<FunctionSpace>;

        bool value(const Vector_t &x, Scalar_t &value) const override { return false; }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            ensure_gradient(g);

            assert(false && "IMPLEMENT ME");
            // if (!assembler_->assemble_vector(x, g)) {
            //     return false;
            // }

            integrate_gradient(x, g);
            return true;
        }

        bool update(const Vector_t &x) override { return true; }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            ensure_hessian(H);

            assert(false && "IMPLEMENT ME");
            // if (!assembler_->assemble_matrix(x, H)) {
            //     return false;
            // }

            integrate_hessian(x, H);
            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            ensure_hessian(H);
            ensure_gradient(g);

            if (!assembler_->assemble(x, H, g)) {
                return false;
            }

            integrate_hessian(x, H);
            integrate_gradient(x, g);
            return true;
        }

        void read(Input &in) override {
            Super::read(in);

            in.get("assembly", *assembler_);

            bool user_defined_mass = false;
            in.get("mass", [this, &user_defined_mass](Input &node) {
                mass_matrix_assembler_->read(node);
                user_defined_mass = true;
            });

            if (!user_defined_mass) {
                auto params = param_list(param("material",
                                               param_list(param("type", "Mass"),                         //
                                                          param("lumped", true),                         //
                                                          param("n_components", this->space()->n_var())  //
                                                          )));
                mass_matrix_assembler_->read(params);
            }

            // if (this->verbose()) {
            //     utopia::out() << "--------------------------------\n";
            //     this->describe(utopia::out().stream());
            //     utopia::out() << "--------------------------------\n";
            // }
        }

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            assembler_->set_environment(env);
            mass_matrix_assembler_->set_environment(env);
        }

        bool init_IVP(Vector_t &x) override {
            this->space()->create_matrix(*mass_matrix_);

            if (empty(x)) {
                this->space()->create_vector(x);
                x.set(0.0);
            }

            this->space()->apply_constraints(x);

            assert(false && "IMPLEMENT ME");
            // if (!mass_matrix_assembler_->assemble_matrix(x, *mass_matrix_)) {
            //     return false;
            // }

            rename("mass_matrix", *mass_matrix_);

            this->space()->create_vector(active_stress_);
            active_stress_.set(0.);

            auto vlo = layout(active_stress_);

            x_old_.zeros(vlo);
            x_older_.zeros(vlo);
            active_stress_old_.zeros(vlo);
            active_stress_older_.zeros(vlo);
            external_force_.zeros(vlo);

            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            const Scalar_t dt = this->delta_time();

            x_older_ = x_old_;
            x_old_ = x;

            active_stress_older_ = active_stress_old_;
            active_stress_old_ = active_stress_;

            active_stress_ = (4. / (dt * dt)) * ((*mass_matrix_) * (2. * x_old_ - x_older_)) - 2. * active_stress_old_ -
                             active_stress_older_ + (4. * external_force_);

            return true;
        }

        NewmarkIntegrator(const std::shared_ptr<FunctionSpace> &space)
            : space_(space),
              assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_(std::make_shared<Matrix_t>()) {}

        inline const std::shared_ptr<FunctionSpace> &space() const { return space_; }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<OmniAssembler_t> assembler_;
        std::shared_ptr<OmniAssembler_t> mass_matrix_assembler_;
        std::shared_ptr<Matrix_t> mass_matrix_;

        Vector_t x_old_, x_older_;
        Vector_t active_stress_, active_stress_old_, active_stress_older_, external_force_;

        void integrate_gradient(const Vector_t &x, Vector_t &g) const {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();
            g -= active_stress_;
            g *= (dt2 / 4.);
            g += ((*mass_matrix_) * x);
            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();
            H *= (dt2 / 4.);
            H += (*mass_matrix_);
            this->space()->apply_constraints(H);
        }

        void ensure_gradient(Vector_t &g) const {
            if (empty(g)) {
                this->space()->create_vector(g);
                rename("g", g);
            } else {
                g *= 0.0;
            }
        }

        void ensure_hessian(Matrix_t &H) const {
            if (empty(H)) {
                this->space()->create_matrix(H);
                rename("H", H);
            } else {
                if (H.is_assembled()) {
                    H *= 0.0;
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP
