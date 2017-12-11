#ifndef UTOPIA_MECHANICS_HPP
#define UTOPIA_MECHANICS_HPP 

#include "utopia.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"

#include <memory>

namespace libMesh {
	class DofMap;
}

namespace utopia {
	class Contact;

	class MechanicsState {
	public:
		void init(const Size &local_size, const Size &global_size);

		DVectord displacement;
		DVectord displacement_increment;
		
		DVectord velocity;
		
		DVectord internal_force;
		DVectord external_force;

		DVectord stress;

		double t;
	};

	class MechanicsContext {
	public:
		//stiffness matrix
		DSMatrixd stiffness_matrix;

		//lumped mass matrix
		DSMatrixd mass_matrix;
		DVectord  inverse_mass_vector;

		void init_mass_matrix(const ProductFunctionSpace<LibMeshFunctionSpace> &V);
	};

	class Friction {
	public:
		Friction()
		: friction_coefficient(0.)
		{}

		//FIXME
		double friction_coefficient;
	};

	class MechIntegrationScheme {
	public:
		virtual ~MechIntegrationScheme() {}

		virtual void apply(
			const double dt,
			const MechanicsContext &mech_ctx,
			const MechanicsState &old,
			MechanicsState &current) = 0;
	};

	class MechWithContactIntegrationScheme {
	public:
		virtual ~MechWithContactIntegrationScheme() {}

		MechWithContactIntegrationScheme(
			const unsigned int dim,
			libMesh::DofMap &dof_map)
		: dim(dim), dof_map(dof_map) 
		{
			linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
		}

		virtual void apply(
			const double dt,
			const MechanicsContext &mech_ctx,
			const Contact  &contact,
			const Friction &friction,
			const MechanicsState &old,
			MechanicsState &current) = 0;

		bool solve(
			const DSMatrixd &K,
			const DVectord &inverse_mass_vector,
			const DVectord &rhs,
			const DVectord &gap,
			const Friction &friction,
			DVectord &sol);

		//FIXME I do not like it
		unsigned int dim;
		libMesh::DofMap &dof_map;
		std::shared_ptr< LinearSolver<DSMatrixd, DVectord> > linear_solver;
	};

	class ImplicitEuler : public MechIntegrationScheme, public MechWithContactIntegrationScheme {
	public:

		ImplicitEuler(
			const unsigned int dim,
			libMesh::DofMap &dof_map)
		: MechWithContactIntegrationScheme(dim, dof_map)
		{ }

		void apply(
			const double dt,
			const MechanicsContext &mech_ctx,
			const MechanicsState &old,
			MechanicsState &current) override;

		void apply(
			const double dt,
			const MechanicsContext &mech_ctx,
			const Contact  &contact,
			const Friction &friction,
			const MechanicsState &old,
			MechanicsState &current) override;
	};

	class ExternalForce {
	public:
		virtual ~ExternalForce() {}
		virtual void eval(const double t, DVectord &result) = 0;
	};

	class ConstantExternalForce : public ExternalForce {
	public:
		inline void eval(const double, DVectord &result) override
		{
			result = value;
		}

		template<class LinearForm>
		void init(const LinearForm &linear_form)
		{
			assemble(linear_form, value);
		}

		DVectord value;
	};

	// class MechanicsSimulation {
	// public:

	// };
}

#endif //UTOPIA_MECHANICS_HPP
