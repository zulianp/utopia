#ifndef UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP
#define UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP 

#include "utopia_SteadyContact.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	class ContactStabilizedNewmark : public ContactSolver<Matrix, Vector> {
	public:
		DEF_UTOPIA_SCALAR(Matrix);

		virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
		{
			if(!this->material().assemble_hessian_and_gradient(x, stiffness_matrix_, internal_force_)) {
				return false;
			}

			auto &O = this->contact().orthogonal_trafo;
			const auto &gap = this->contact().gap;

			//predictor step
			pred_ = O * utopia::min(O * (dt_ * velocity_), gap);

			gradient = internal_mass_matrix_ * pred_ + (dt_* dt_/2.) * (2. * external_force_ - internal_force_);
			hessian  = internal_mass_matrix_ + (dt_ * dt_/4.) * stiffness_matrix_;

			Fcon_m_F_ = (-2./dt_) * (internal_mass_matrix_ * (pred_ - (x - old_x_)));
			apply_zero_boundary_conditions(this->space()[0].dof_map(), Fcon_m_F_);
			velocity_ += e_mul(inverse_mass_vector_, Fcon_m_F_);

			return true;
		}

		virtual void finalize() override
		{
		
		}

	private:
		Scalar dt_;
		Matrix stiffness_matrix_;
		Matrix internal_mass_matrix_;
		Vector internal_force_;
		Vector external_force_;
		Vector velocity_;
		Vector inverse_mass_vector_;
		Vector old_x_;

		//buffers
		Vector pred_;
		Vector Fcon_m_F_;

	};

}


#endif //UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP