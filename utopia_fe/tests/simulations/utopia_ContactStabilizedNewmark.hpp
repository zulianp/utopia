// #ifndef UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP
// #define UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP 

// namespace utopia {


// 	template<class Matrix, class Vector, class Eqs>
// 	class ContactStabilizedNewmark : public Function<Matrix, Vector> {
// 	public:
// 		DEF_UTOPIA_SCALAR(Matrix)

// 		ContactStabilizedNewmark()
// 		: eqs_(eqs), first_(true)
// 		{}

// 		virtual bool value(const Vector &, Scalar &) const override
// 		{
// 			return false;
// 		}

// 		virtual bool gradient(const Vector &x, Vector &result) const override
// 		{
// 			result = g_;
// 			return true;
// 		}

// 		virtual bool hessian(const Vector &, Matrix &result) const override
// 		{
// 			result = H_;
// 			return true;
// 		}

// 		void update_contact(const Vector &x)
// 		{

// 		}

// 		//compute sub-differential
// 		virtual bool update(const Vector &x) override { 
// 			update_contact(x);

// 			const auto &T   = contact_.complete_transformation;
// 			const auto &O   = contact_.orthogonal_trafo;
// 			const auto &gap = contact_.gap;

// 			//predictor step
// 			DVectord pred = O * utopia::min(O * (dt * velocity), gap);
// 			DVectord rhs  = internal_mass_matrix * pred + (dt*dt/2.) * (2. * external_force - internal_force);
// 			DSMatrixd K   = internal_mass_matrix + (dt*dt/4.) * stiffness_matrix;

// 			DVectord sol_c = local_zeros(local_size(external_force));
// 			DVectord rhs_c = transpose(T) * rhs;
// 			DSMatrixd K_c  = transpose(T) * K * T;

// 			auto u = fe_function(*spaces[0]);
// 			apply_boundary_conditions(u, K_c, rhs_c);

// 			if(has_friction) {
// 				solve_with_friction(K_c, rhs_c, sol_c);
// 			} else {
// 				solve_without_friction(K_c, rhs_c, sol_c);
// 			}

			
// 			return true;
// 		}

// 		void time_step_begin()
// 		{

// 		}

// 		void time_step_end()
// 		{
// 			// displacement_increment = T * sol_c;		
// 			// total_displacement += displacement_increment;


// 			// new_internal_force = stiffness_matrix * total_displacement;
// 			// apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);

// 			// DVectord Fcon_m_F = (-2./dt) * (internal_mass_matrix * (pred - displacement_increment));

// 			// apply_zero_boundary_conditions(spaces[0]->dof_map(), Fcon_m_F);
// 			// DVectord vel_inc  = e_mul(inverse_mass_vector, Fcon_m_F);
// 			// velocity += vel_inc;
// 		}

// 		void reset()
// 		{
// 			first_ = true;
// 		}

// 	private:
// 		bool first_;
// 		Matrix H_;
// 		Vector g_;

// 		Contact contact_;

// 		//states
// 		MechanicsState current_state;
// 		MechanicsState previous_state;

// 		MechanicsContext mech_ctx;

// 		//additional vectors
// 		DVectord Fcon_m_F_;
// 		DSMatrixd internal_mass_matrix_;
// 	};

// }


// #endif //UTOPIA_CONTACT_STABILIZED_NEWMARK_HPP