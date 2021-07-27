#ifndef UTOPIA_NEWMARK_HPP
#define UTOPIA_NEWMARK_HPP

#include "utopia_ContactSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class Newmark : public ContactSolver<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef typename ContactSolver<Matrix, Vector>::FunctionSpaceT FunctionSpaceT;

        Newmark(const std::shared_ptr<FunctionSpaceT> &V,
                const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
                const Scalar dt,
                const ContactParams &params)
            : ContactSolver<Matrix, Vector>(V, material, params), dt_(dt) {}

        // virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            if (!this->material().assemble_hessian_and_gradient(x, stiffness_matrix_, internal_force_)) {
                return false;
            }

            // hessian  = (4./(dt_*dt_)) * internal_mass_matrix_ + stiffness_matrix_;
            // gradient = internal_force_ + 4./(dt_*dt_) * (internal_mass_matrix_ * x) - forcing_term_;

            auto dt2 = dt_ * dt_;
            hessian = internal_mass_matrix_ + (dt2 / 4.) * stiffness_matrix_;
            gradient = (dt2 / 4.) * (internal_force_ + 4. / (dt2) * (internal_mass_matrix_ * x) - forcing_term_);
            return true;
        }

        void initialize() override { ContactSolver<Matrix, Vector>::initialize(); }

        void initial_condition() {
            auto &V = this->space();
            auto u = trial(V);
            auto v = test(V);

            utopia::assemble(inner(u, v) * dX, internal_mass_matrix_);
            Vector mass_vector = sum(internal_mass_matrix_, 1);
            internal_mass_matrix_ = diag(mass_vector);

            auto n_local = local_size(internal_mass_matrix_).get(0);
            external_force_ = local_zeros(n_local);
            internal_force_old_ = local_zeros(n_local);
            internal_force_older_ = local_zeros(n_local);

            x_old_ = local_zeros(n_local);
            x_older_ = local_zeros(n_local);
            forcing_term_ = local_zeros(n_local);

            t_ = 0.;

            if (this->external_force_fun()) {
                this->external_force_fun()->eval(t_, external_force_);
            }
        }

        void next_step() override {
            x_older_ = x_old_;
            x_old_ = this->displacement();

            internal_force_older_ = internal_force_old_;
            internal_force_old_ = internal_force_;

            forcing_term_ = (4. / (dt_ * dt_)) * (internal_mass_matrix_ * (2. * x_old_ - x_older_)) -
                            2. * internal_force_old_ - internal_force_older_ + (4. * external_force_);
        }

        void finalize() override {}

    private:
        Scalar dt_;
        Scalar t_;

        // operators
        Matrix stiffness_matrix_;
        Matrix internal_mass_matrix_;

        // old displacement
        Vector x_old_;
        Vector x_older_;

        // forces
        Vector internal_force_;
        Vector internal_force_old_;
        Vector internal_force_older_;

        Vector external_force_;
        Vector forcing_term_;
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_HPP