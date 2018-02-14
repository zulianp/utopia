#ifndef UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP
#define UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP 

#include "utopia_SteadyContact.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	class ContactStabilizedNewmark : public ContactSolver<Matrix, Vector> {
	public:
		DEF_UTOPIA_SCALAR(Matrix);
		typedef typename ContactSolver<Matrix, Vector>::FunctionSpaceT FunctionSpaceT;

		ContactStabilizedNewmark(
			const std::shared_ptr<FunctionSpaceT> &V,
			const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
			const Scalar dt,
			const ContactParams &params
		)
		: ContactSolver<Matrix, Vector>(V, material, params), dt_(dt), is_new_time_step_(true)
		{}

		// virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
		virtual bool assemble_hessian_and_gradient(Vector &x, Matrix &hessian, Vector &gradient) override
		{
			if(!this->material().assemble_hessian_and_gradient(x, stiffness_matrix_, internal_force_)) {
				return false;
			}

			if(is_new_time_step_) {
				auto &O = this->contact().orthogonal_trafo;
				pred_ = O * utopia::min(O * (dt_ * velocity_old_), this->contact().gap);
				is_new_time_step_ = false;
			}

			hessian  = internal_mass_matrix_ + ((dt_*dt_)/4.) * stiffness_matrix_;
			gradient = ((dt_*dt_)/4.) * internal_force_ + (internal_mass_matrix_ * (x - pred_)) - forcing_term_;
			return true;
		}

		virtual void initialize() override
		{
			ContactSolver<Matrix, Vector>::initialize();
		}

		void initial_condition()
		{
			auto &V = this->space();
			auto u = trial(V);
			auto v = test(V);

			utopia::assemble(inner(u, v) * dX, internal_mass_matrix_);
			Vector mass_vector = sum(internal_mass_matrix_, 1);
			internal_mass_matrix_ = diag(mass_vector);
			inverse_mass_vector_ = 1./mass_vector;

			auto n_local 		= local_size(internal_mass_matrix_).get(0);
			external_force_		= local_zeros(n_local);
			internal_force_old_ = local_zeros(n_local);

			x_old_        = local_zeros(n_local);
			forcing_term_ = local_zeros(n_local);
			velocity_     = local_zeros(n_local);
			velocity_old_ = local_zeros(n_local);

			t_ = 0.;

			if(this->external_force_fun()) {
				this->external_force_fun()->eval(t_, external_force_);
			}
		}

		inline const Vector &velocity() const
		{
			return velocity_;
		}

		void update_velocity()
		{
			velocity_inc_ = (-2./dt_) * (internal_mass_matrix_ * (x_old_ -  this->displacement() + pred_));
			apply_zero_boundary_conditions(this->space()[0].dof_map(), velocity_inc_);			
			velocity_ = velocity_old_ + e_mul(inverse_mass_vector_, velocity_inc_);
		}

		void next_step() override
		{	
			update_velocity();

			x_old_ = this->displacement();
			internal_force_old_ = internal_force_;
			velocity_old_ = velocity_;

			//external_force_ = (external_force_old + external_force_current)/2
			forcing_term_ = internal_mass_matrix_ * x_old_ + (dt_*dt_/4.) * (2. * external_force_ - internal_force_old_);
			is_new_time_step_ = true;
		}

		virtual void finalize() override
		{
		
		}

	private:
		Scalar dt_;
		Scalar t_;
		
		//operators
		Matrix stiffness_matrix_;
		Matrix internal_mass_matrix_;
		Vector inverse_mass_vector_;
	
		//old displacement
		Vector x_old_;

		Vector velocity_;
		Vector velocity_old_;


		//forces
		Vector internal_force_;	
		Vector internal_force_old_;

		Vector external_force_;
		Vector forcing_term_;
		Vector velocity_inc_;

		Vector pred_;
		bool is_new_time_step_;

	};

}


#endif //UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP