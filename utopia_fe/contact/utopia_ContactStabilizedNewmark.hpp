#ifndef UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP
#define UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP

#include "utopia_fe_base.hpp"
#include "utopia_libmesh_Types.hpp"

#ifndef WITH_TRILINOS_ALGEBRA

#include "utopia_ContactSolver.hpp"
#include "utopia_ContactSystem.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class ContactStabilizedNewmark : public ContactSolver<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        using super = ContactSolver<Matrix, Vector>;
        using FunctionSpaceT = typename super::FunctionSpaceT;
        using ContactT = typename super::ContactT;

        ContactStabilizedNewmark(const std::shared_ptr<FunctionSpaceT> &V,
                                 const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
                                 const Scalar dt,
                                 const ContactParams &params)
            : ContactSolver<Matrix, Vector>(V, material, params), dt_(dt), density_(1.), is_new_time_step_(true) {
            // c_sys_ = std::make_shared<ContactSystem>(V->subspace(0).equation_systems_ptr(),
            // V->subspace(0).equation_system().number());
        }

        // virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            if (!this->material().assemble_hessian_and_gradient(x, stiffness_matrix_, internal_force_)) {
                return false;
            }

            if (is_new_time_step_) {
                // copy for time-adaptivity for changing the predictor for different dt

                O_copy_ = this->contact().orthogonal_trafo();
                // O_copy_   = this->contact().complete_transformation();
                gap_copy_ = this->contact().gap();
                update_predictor();

                is_new_time_step_ = false;
            }

            hessian = internal_mass_matrix_ + ((dt_ * dt_ * density_) / 4.) * stiffness_matrix_;
            gradient =
                ((dt_ * dt_ * density_) / 4.) * internal_force_ + (internal_mass_matrix_ * (x - pred_)) - forcing_term_;
            return true;
        }

        bool stress(const Vector &x, Vector &result) override {
            result =
                ((dt_ * dt_ * density_) / 4.) * internal_force_ + (internal_mass_matrix_ * (x - pred_)) - forcing_term_;
            return true;
        }

        void initialize() override { ContactSolver<Matrix, Vector>::initialize(); }

        void initial_condition(const Scalar density, const UVector &initial_velocity = UVector()) {
            auto &V = this->space();
            auto u = trial(V);
            auto v = test(V);

            utopia::assemble(inner(u, v) * dX, internal_mass_matrix_);
            // internal_mass_matrix_ *= density;
            initial_condition(internal_mass_matrix_, initial_velocity);
            density_ = density;
        }

        void initial_condition(const USparseMatrix &mass_matrix, const UVector &initial_velocity = UVector()) {
            internal_mass_matrix_ = mass_matrix;
            Vector mass_vector = sum(internal_mass_matrix_, 1);
            internal_mass_matrix_ = diag(mass_vector);
            inverse_mass_vector_ = 1. / mass_vector;

            auto n_local = local_size(internal_mass_matrix_).get(0);
            external_force_ = local_zeros(n_local);
            internal_force_old_ = local_zeros(n_local);

            x_old_ = local_zeros(n_local);
            forcing_term_ = local_zeros(n_local);
            velocity_ = local_zeros(n_local);

            if (empty(initial_velocity)) {
                velocity_old_ = local_zeros(n_local);
            } else {
                velocity_old_ = initial_velocity;
            }

            velocity_inc_ = local_zeros(n_local);
            pred_ = local_zeros(n_local);
            gap_copy_ = local_zeros(n_local);

            t_ = 0.;

            if (this->external_force_fun()) {
                this->external_force_fun()->eval(t_, external_force_);
            }
        }

        void initial_condition(const UVector &initial_velocity = UVector()) {
            auto &V = this->space();
            auto u = trial(V);
            auto v = test(V);

            utopia::assemble(inner(u, v) * dX, internal_mass_matrix_);
            initial_condition(internal_mass_matrix_, initial_velocity);
        }

        inline const Vector &velocity() const { return velocity_; }

        void update_predictor() {
            pred_ = O_copy_ * utopia::min(dt_ * (transpose(O_copy_) * velocity_old_), gap_copy_);
        }

        void update_velocity() {
            velocity_inc_ = (-2. / dt_) * (internal_mass_matrix_ * (x_old_ - this->displacement() + pred_));
            apply_zero_boundary_conditions(this->space()[0].dof_map(), velocity_inc_);
            velocity_ = velocity_old_ + e_mul(inverse_mass_vector_, velocity_inc_);
        }

        void update_forcing_term() {
            // if(this->external_force_fun()) {
            // 	this->external_force_fun()->eval(t_, external_force_);
            // }

            // FIXME is density_ * external_force wrong?
            forcing_term_ = internal_mass_matrix_ * x_old_ +
                            ((dt_ * dt_ * density_) / 4.) * (2. * external_force_ - internal_force_old_);
        }

        void next_step() override {
            update_velocity();
            update_contact_system();

            x_old_ = this->displacement();
            internal_force_old_ = internal_force_;
            velocity_old_ = velocity_;

            // external_force_ = (external_force_old + external_force_current)/2

            update_forcing_term();
            is_new_time_step_ = true;
            t_ += dt_;
        }

        void update_contact_system() {
            if (!c_sys_) return;

            MechanicsState state;
            state.displacement = this->displacement();
            state.displacement_increment = this->displacement() - x_old_;

            state.velocity = velocity_;

            state.internal_force = internal_force_;
            state.external_force = external_force_;

            state.stress = external_force_ - internal_force_;
            c_sys_->update(state, this->contact(), dt_);
        }

        void finalize() override {}

        void set_dt_and_update(const Scalar dt) {
            if (approxeq(dt, dt_)) {
                return;
            }

            dt_ = dt;

            update_forcing_term();
            update_predictor();
        }

    private:
        Scalar dt_;
        Scalar t_;
        Scalar density_;

        // operators
        Matrix stiffness_matrix_;
        Matrix internal_mass_matrix_;
        Vector inverse_mass_vector_;

        // old displacement
        Vector x_old_;

        Vector velocity_;
        Vector velocity_old_;

        // forces
        Vector internal_force_;
        Vector internal_force_old_;

        Vector external_force_;
        Vector forcing_term_;
        Vector velocity_inc_;

        Vector pred_;
        Vector gap_copy_;
        Matrix O_copy_;
        bool is_new_time_step_;

        std::shared_ptr<ContactSystem> c_sys_;
    };

}  // namespace utopia

#endif  // WITH_TRILINOS_ALGEBRA
#endif  // UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP
