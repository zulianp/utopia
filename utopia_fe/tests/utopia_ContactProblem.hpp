#ifndef UTOPIA_CONTACT_PROBLEM_HPP
#define UTOPIA_CONTACT_PROBLEM_HPP 

#include "utopia_fe_core.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "moonolith_communicator.hpp"
#include "utopia.hpp"

namespace utopia {
	class ContactProblem {
	public:
		typedef utopia::FESpace<LibMeshTraits<libMesh::Real> > FESpaceT;
		typedef std::shared_ptr<FESpaceT> FESpacePtr;
		typedef utopia::LibMeshFEContext<libMesh::LinearImplicitSystem> FEContextT;
		typedef std::shared_ptr<FEContextT> FEContextPtr;

		//discretization
		FEContextPtr context_ptr;
		std::shared_ptr<libMesh::MeshBase> mesh;
		// std::shared_ptr<libMesh::MeshBase> displaced_mesh;
		std::vector<FESpacePtr> spaces;

		//material quantities
		DSMatrixd stiffness_matrix;
		DSMatrixd neumann_matrix;
		DSMatrixd mass_matrix;
		DVectord force;
		DVectord displacement;
		DVectord old_displacement;
		DVectord total_displacement;

		DVectord velocity;
		DVectord old_velocity;

		DVectord acceleration;
		DVectord old_acceleration;

		DVectord internal_force;

		//contact quantities
		DSMatrixd coupling;
		DSMatrixd transfer_operator;
		DSMatrixd orthogonal_trafo;
		DVectord normals;
		DVectord is_contact_node;
		DVectord weighted_gap;
		DVectord gap;
		DSMatrixd boundary_mass_inv;
		DVectord normal_stress;

		double search_radius;
		std::vector< std::pair<int, int> > contact_pair_tags;

		moonolith::Communicator comm;
		std::shared_ptr< LinearSolver<DSMatrixd, DVectord> > linear_solver;

		int iteration;
		bool verbose;

		class ElasticityBoundaryConditions {
		public:
			virtual ~ElasticityBoundaryConditions() {}
			virtual void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy) = 0;
			virtual void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) = 0;
		};

		std::shared_ptr< ElasticityBoundaryConditions > bc_ptr;

		void step(const double dt = 1.0);

		void init(
			const libMesh::LibMeshInit &init, 
			const std::shared_ptr<libMesh::MeshBase> &mesh,
			const std::shared_ptr< ElasticityBoundaryConditions > &bc_ptr,
			std::vector< std::pair<int, int> > contact_pair_tags,
			double search_radius
		);

		void save(const std::string &output_dir = ".");

		ContactProblem();
	private:

		void init_discretization();
		void init_material();
		void compute_contact_conditions();
		void init_material_2d();
		void init_material_3d();
		void apply_displacement(const DVectord &displacement);
	};
}

#endif //UTOPIA_CONTACT_PROBLEM_HPP
