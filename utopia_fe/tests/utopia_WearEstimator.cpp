#include "utopia_WearEstimator.hpp"

#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"
#include "utopia_Socket.hpp"
#include "utopia_ContactProblem.hpp"

#include "libmesh/exodusII_io.h"
#include <algorithm>


#include <memory>
#include <unordered_set>


namespace utopia {



	class SubdomainAffineTransform {
	public:
		libMesh::TensorValue<double> linear_transformation;
		libMesh::VectorValue<double> translation;
		int subdomain_id;

		SubdomainAffineTransform()
		: subdomain_id(0)
		{}

		void apply(libMesh::MeshBase &mesh) const
		{
			// for(auto n_it = mesh.local_nodes_begin(); n_it != mesh.local_nodes_end(); ++n_it) { 
			// 	auto &node = **n_it;
			// 	node = linear_transformation * node;
			// 	node += translation;
			// }


			std::unordered_set<libMesh::dof_id_type> transformed;

			auto m_begin = mesh.active_local_elements_begin();
			auto m_end   = mesh.active_local_elements_end();

			for(auto m_it = m_begin; m_it != m_end; ++m_it) { 
				auto &e = **m_it;

				if(e.subdomain_id() != subdomain_id) continue;

				for(int i = 0; i < e.n_nodes(); ++i) {
					auto &node = e.node_ref(i);
					auto ret = transformed.insert(node.id());

					if(ret.second) {
						node = linear_transformation * node;
						node += translation;
					}
				}
			}
		}

		void make_rotation(const int dim, const double angle, const char axis)
		{
			translation.zero();

			switch(dim) {
				case 2:
				{	
					assert(axis == 'y' || axis == 'Y');
					make_rotation_2(angle, linear_transformation);
					return;	
				}
				case 3:
				{
					make_rotation_3(angle, axis, linear_transformation);
					return;
				}
				default: 
				{
					assert(false);
					return;
				}
			}
		}

		static void make_rotation_3(const double angle, const char axis, libMesh::TensorValue<double> &result)  
		{

			result(0, 0) = result(1, 1) = result(2, 2) = 1.;

			if ((axis == 'x') || (axis == 'X')) {
				result(1, 1) =  cos(angle);   
				result(1, 2) = -sin(angle);   			
				result(2, 1) =  sin(angle);   		
				result(2, 2) =  cos(angle);   			
			} else
				// set rotation around y axis:
				if ((axis == 'y') || (axis == 'Y')) {
					result(0, 0) =  cos(angle);   
					result(0, 2) =  sin(angle);   
					result(2, 0) = -sin(angle);   
					result(2, 2) =  cos(angle);   
				} else
					// set rotation around z axis:
					if ((axis == 'z') || (axis == 'Z')) {
						result(0, 0) =  cos(angle);   
						result(0, 1) = -sin(angle);   
						result(1, 0) =  sin(angle);   
						result(1, 1) =  cos(angle);   
					}
		}

		static void make_rotation_2(const double angle, libMesh::TensorValue<double> &result)
		{
			result(0, 0) = cos(angle);
			result(0, 1) = -sin(angle);
			result(1, 0) = sin(angle);
			result(1, 1) = cos(angle);
		}

	};

	class SubdomainCompositeTransform {
	public:
		void apply(libMesh::MeshBase &mesh) const
		{
			for(const auto &t : sub_transforms) {
				t.apply(mesh);
			}
		}

		std::vector<SubdomainAffineTransform> sub_transforms;
	};

	class GaitCycle {
	public:

		class SetUp : 
			public ContactProblem::ElasticityBoundaryConditions,
			public ContactProblem::ElasticityForcingFunction {

		public:

			SetUp() {
				contact_flags = {{2, 1}};
				search_radius = 0.1;
				dt = .01;
				
				params.set_mu(1, 2.0);
				params.set_lambda(1, 2.0);

				params.set_mu(2, 5.0);
				params.set_lambda(2, 5.0);

				wear_coefficient = 7e-3; //mild-steel
			}

			void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy) override 
			{
				strong_enforce( boundary_conditions(ux == coeff(0.), {3}) );
				strong_enforce( boundary_conditions(uy == coeff(0.), {3}) );

				strong_enforce( boundary_conditions(ux == coeff(0.),  {4}) );
				strong_enforce( boundary_conditions(uy == coeff(0.1), {4}) );
			}

			void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) override
			{
				strong_enforce( boundary_conditions(ux == coeff(0.), {3}) );
				strong_enforce( boundary_conditions(uy == coeff(0.), {3}) );
				strong_enforce( boundary_conditions(uz == coeff(0.), {3}) );

				strong_enforce( boundary_conditions(ux == coeff(0.), {4}) );
				strong_enforce( boundary_conditions(uz == coeff(0.), {4}) );
			}

			int block_id() const override
			{
				return 2;
			}

			virtual void fill(libMesh::DenseVector<libMesh::Real> &v) override
			{
				v.zero();
				// v(1) = 0.1;
			}

			LameeParameters params;
			std::vector<std::pair<int, int> > contact_flags;
			double search_radius;
			double dt;
			double wear_coefficient;
		};

		std::shared_ptr<libMesh::MeshBase> mesh;
		SubdomainAffineTransform transform;
		DSMatrixd mass_matrix;
		DSMatrixd stiffness_matrix;
		ContactProblem contact_problem;

		void apply(libMesh::LibMeshInit &init, const std::size_t n_configurations)
		{
			auto set_up = std::make_shared<SetUp>();
			ContactProblem p;
			p.has_friction = false;

			for(std::size_t i = 0; i < n_configurations; ++i) {
				transform.apply(*mesh);
				plot_mesh(*mesh, "mesh");
				p.params = set_up->params;
				p.must_apply_displacement = false;
				p.init(init, mesh, set_up, set_up, set_up->contact_flags, set_up->search_radius);
				p.step(set_up->dt);
				p.save(set_up->dt, i, "wear_simulation" + std::to_string(i) + ".e");	
			}
		}
	};

	class WearEstimator {
	public:

		void run(libMesh::LibMeshInit &init, const std::shared_ptr<libMesh::MeshBase> &mesh) {
			GaitCycle gait_cycle;
			gait_cycle.mesh = mesh;

			const int dim = mesh->mesh_dimension();

			const int n_configurations = 2;
			const double d_angle = 60. * (M_PI/180.) / n_configurations;

			auto &t = gait_cycle.transform;
			t.make_rotation(dim, d_angle, 'y');
			t.subdomain_id = 1;

			// m1 = LameeParamters::convert_to_lamee_params(220., 0.31);
			// m2 = LameeParamters::convert_to_lamee_params(1.1, 0.42);

			gait_cycle.apply(init, n_configurations);
		}
	};

	void run_wear_test(libMesh::LibMeshInit &init)
	{

		std::cout << "[run_wear_test]" << std::endl;

		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	
		mesh->read("../data/wear_2.e");

		WearEstimator we;
		we.run(init, mesh);
	}
}
