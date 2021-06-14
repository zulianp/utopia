#ifndef UTOPIA_FSI_HESCH_2014_HPP
#define UTOPIA_FSI_HESCH_2014_HPP

#include "utopia_FEModelFunction.hpp"

namespace utopia {

    /**
     * @brief from https://www.sciencedirect.com/science/article/abs/pii/S0045782514001893?via%3Dihub
     *
     */
    template <class FunctionSpace>
    class FSIHesch2014 final : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Environment_t = utopia::Environment<FunctionSpace>;

        using FunctionPtr_t = std::shared_ptr<FEFunctionInterface<FunctionSpace>>;
        using Field_t = utopia::Field<FunctionSpace>;

        void create_solution_vector(Vector_t &fluid_x) override { fluid_->create_solution_vector(fluid_x); }

        void apply_constraints(Vector_t &fluid_x) const override { fluid_->apply_constraints(fluid_x); }
        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            fluid_->set_environment(env);
            solid_->set_environment(env);
        }

        const std::shared_ptr<Matrix_t> &mass_matrix() const override { return fluid_->mass_matrix(); }

        bool assemble_mass_matrix() override { return fluid_->assemble_mass_matrix(); }
        bool assemble_mass_matrix(Matrix_t &mass_matrix) override { return fluid_->assemble_mass_matrix(mass_matrix); }

        const std::shared_ptr<FunctionSpace> &space() const override { return fluid_->space(); }

        bool is_time_dependent() const override { return fluid_->is_time_dependent(); }
        bool is_linear() const override { return fluid_->is_linear() && solid_->is_linear(); }

        void must_apply_constraints_to_assembled(const bool val) override {
            fluid_->must_apply_constraints_to_assembled(val);
        }

        bool report_solution(const Vector_t &fluid_x) override {
            return fluid_->report_solution(fluid_x) && solid_->report_solution(displacement_->data());
        }

        void update_FSI_quantities(const Vector_t &fluid_x) const {
            Vector_t solid_velocity;
            transfer_.apply(fluid_x, *diff_velocity_);

            solid_->time_derivative(displacement_->data(), solid_velocity);

            *diff_velocity_ -= solid_velocity;
        }

        bool update_IVP(const Vector_t &fluid_x) override { return fluid_->update_IVP(fluid_x); }
        bool setup_IVP(Vector_t &fluid_x) override { return fluid_->setup_IVP(fluid_x); }
        bool is_IVP_solved() override { return fluid_->is_IVP_solved(); }

        void update_transfer() {
            solid_->displace(displacement_->data());
            transfer_.init(*fluid_->space(), *solid_->space());
            solid_->displace(-displacement_->data());
        }

        void setup_transfer() {
            FETransferOptions opts;
            opts.n_var = solid_->n_var();
            transfer_.set_options(opts);
        }

        FSIHesch2014(const FunctionPtr_t &fluid, const FunctionPtr_t &solid)
            : fluid_(fluid), solid_(solid), displacement_(std::make_shared<Field_t>("displacement", solid_->space())) {
            displacement_->zero();

            H_solid_ = utopia::make_unique<Matrix_t>();
            g_solid_ = utopia::make_unique<Vector_t>();
            diff_velocity_ = utopia::make_unique<Vector_t>();
        }

        bool value(const Vector_t & /*point*/, Scalar_t &value) const override {
            value = 123456789;
            return false;
        }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(g);
            assert(false);
            return false;
        }

        bool update(const Vector_t &x) override {
            UTOPIA_UNUSED(x);
            return true;
        }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(H);
            assert(false);
            return false;
        }

        bool hessian_and_gradient(const Vector_t &x_FSI, Matrix_t &H_FSI, Vector_t &g_FSI) const override {
            ensure_hessian(H_FSI);
            ensure_gradient(g_FSI);

            if (!fluid_->hessian_and_gradient(x_FSI, H_FSI, g_FSI)) {
                return false;
            }

            update_FSI_quantities(x_FSI);

            if (!solid_->hessian_and_gradient(displacement_->data(), *H_solid_, *g_solid_)) {
                return false;
            }

            Matrix_t H_temp;
            Vector_t g_temp;
            g_temp -= (*H_solid_) * (*diff_velocity_);

            transfer_.apply_transpose(*g_solid_, g_temp);
            transfer_.apply(*H_solid_, H_temp);

            H_FSI += H_temp;
            g_FSI += g_temp;

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H_FSI);
                this->space()->apply_zero_constraints(g_FSI);
            }

            return true;
        }

    private:
        FunctionPtr_t fluid_;
        FunctionPtr_t solid_;
        std::shared_ptr<Field_t> displacement_;
        std::unique_ptr<Matrix_t> H_solid_;
        std::unique_ptr<Vector_t> g_solid_;
        std::unique_ptr<Vector_t> diff_velocity_;

        FETransfer<FunctionSpace> transfer_;

        bool must_apply_constraints_{true};

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

#endif  // UTOPIA_FSI_HESCH_2014_HPP
