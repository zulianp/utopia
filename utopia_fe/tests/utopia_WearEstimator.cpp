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

#include "libmesh/exodusII_io.h"
#include <algorithm>


#include <memory>
#include <unordered_set>


namespace utopia {

	class ModelParamters {
	public:	
		//lamee paramters
		double mu, lambda;

		static ModelParamters convert_to_lamee_params(const double young_modulus, double poisson_ratio)
		{	
			ModelParamters ret;
			return ret;
		}
	};

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

		std::vector<std::shared_ptr<SubdomainCompositeTransform>> transform;
		DSMatrixd mass_matrix;
		DSMatrixd stiffness_matrix;
	};


	class WearEstimator {
	public:


	};

	void run_wear_test(libMesh::LibMeshInit &init)
	{

		std::cout << "[run_wear_test]" << std::endl;

		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	
		mesh->read("../data/wear_2.e");

		// const unsigned int n = 10;
		// libMesh::MeshTools::Generation::build_square(*mesh,
		// 	n, n,
		// 	0, 1,
		// 	0, 1.,
		// 	libMesh::TRI3);

		const int dim = mesh->mesh_dimension();

		GaitCycle gc;
		gc.transform.push_back(std::make_shared<SubdomainCompositeTransform>());
		auto &t = *gc.transform.back();



		t.sub_transforms.emplace_back();
		t.sub_transforms.back().make_rotation(dim, 0.2, 'y');
		t.sub_transforms.back().subdomain_id = 1;
		t.apply(*mesh);


		ModelParamters m1, m2;
		m1 = ModelParamters::convert_to_lamee_params(220., 0.31);
		m2 = ModelParamters::convert_to_lamee_params(1.1, 0.42);

		plot_mesh(*mesh, "mesh");

	}
}
