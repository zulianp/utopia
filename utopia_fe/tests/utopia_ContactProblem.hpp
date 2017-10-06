#ifndef UTOPIA_CONTACT_PROBLEM_HPP
#define UTOPIA_CONTACT_PROBLEM_HPP 

#include "utopia_fe_core.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "moonolith_communicator.hpp"
#include "utopia.hpp"

#include "libmesh/nemesis_io.h"

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
		DSMatrixd internal_mass_matrix;
		DSMatrixd constrained_mass_matrix;
		DVectord external_force;
		DVectord displacement_increment;
		DVectord old_displacement_increment;
		DVectord total_displacement;
		DVectord rays;

		DVectord velocity;
		// DVectord old_velocity;

		// DVectord acceleration;
		// DVectord old_acceleration;

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
		DVectord new_internal_force;

		double search_radius;
		std::vector< std::pair<int, int> > contact_pair_tags;

		moonolith::Communicator comm;
		std::shared_ptr< LinearSolver<DSMatrixd, DVectord> > linear_solver;

		int iteration;
		bool verbose;
		bool dynamic_contact;


		std::vector<int> var_num_aux;



		class ElasticityBoundaryConditions {
		public:
			virtual ~ElasticityBoundaryConditions() {}
			virtual void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy) = 0;
			virtual void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) = 0;
		};

		class ElasticityForcingFunction {
		public:
			virtual ~ElasticityForcingFunction() {}
			virtual int block_id() const
			{
				return -1;
			}

			virtual void fill(libMesh::DenseVector<libMesh::Real> &v) 
			{
				v.zero();
			}
		};

		std::shared_ptr<ElasticityBoundaryConditions> bc_ptr;
		std::shared_ptr<ElasticityForcingFunction> ff_ptr;
		std::shared_ptr<libMesh::Nemesis_IO> output;
		// std::shared_ptr<libMesh::ExodusII_IO> output;

		typedef struct {
			double t;
			double kinetic_energy;
			double elastic_energy;
			double potential_energy;
			double contact_energy;
		} Energy;

		std::vector<Energy> energy;

		void compute_energy(const double dt);
		void save_energy(const std::string &path);


		void step(const double dt = 1.0);

		void init(
			const libMesh::LibMeshInit &init, 
			const std::shared_ptr<libMesh::MeshBase> &mesh,
			const std::shared_ptr< ElasticityBoundaryConditions > &bc_ptr,
			const std::shared_ptr< ElasticityForcingFunction > &ff_ptr,
			std::vector< std::pair<int, int> > contact_pair_tags,
			double search_radius
		);

		void save(const double dt = 1.0, const std::string &output_dir = ".");
		inline void set_dynamic_contact(const bool val)
		{
			dynamic_contact = val;
		}

		void is_inpulse(const bool val)
		{
			is_inpulse_ = val;
		}

		ContactProblem();
	private:
		void compute_normal_stress(const double dt);
		void init_discretization();
		void init_material();
		void compute_contact_conditions();
		void init_material_2d();
		void init_material_3d();
		void apply_displacement(const DVectord &displacement);
		void init_aux_system();
		void contact_stabilized_newmark(const double dt);
		void implicity_euler(const double dt);
		void classic_newmark(const double dt);
		void classic_newmark_with_contact(const double dt);
		void classic_newmark_beta(const double dt);
		void contact_stabilized_newmark_monolithic(const double dt);

		bool is_inpulse_;
	};
}

#endif //UTOPIA_CONTACT_PROBLEM_HPP
