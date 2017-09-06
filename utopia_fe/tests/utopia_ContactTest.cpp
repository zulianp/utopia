#include "utopia_ContactTest.hpp"
#include "utopia_assemble_contact.hpp"


#include "utopia.hpp"

//fe extension
#include "utopia_fe_core.hpp"
#include "MortarAssembler.hpp"
#include "ParMortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>

#include "utopia_ContactSimParams.hpp"
#include "libmesh/linear_partitioner.h"
#include "LibmeshContactForMoose.hpp"
#include "LibmeshTransferForMoose.hpp"
#include "LibmeshTransferForMooseReverse.hpp"

#include "utopia_Polygon.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include "moonolith_profiler.hpp"
#include "utopia_ContactProblem.hpp"

#include <iostream>

using namespace libMesh;
using std::make_shared;
using std::shared_ptr;

namespace utopia {


	class ExampleProblem3D : public ContactProblem::ElasticityBoundaryConditions {
	public:
		ExampleProblem3D() {
			mesh_file = "../data/multibody.e";
			top_boundary_tag = 1;
			bottom_boundar_tag = 2;
			contact_flags = {{12, 13}, {14, 13}};
			is_multibody = true;
			search_radius = 0.05;

			//---------------------------------------------------
			// mesh->read("../data/hertz_530.e");
			// mesh->read("../data/quasi_signorini_4593.e");
			// mesh->read("../data/quasi_signorini_fine.e");
			// mesh->read("../data/quasi_signorini_fine_surface_both.e");
			// mesh->read("../data/quasi_signorini_ultra_fine_surface_both.e");
			// mesh->read("../data/two_rocks_26653.e");
			// mesh->read("../data/quasi_signorini_526.e");
			// mesh->read("../data/quasi_signorini_248322.e");
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &)  override {}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) override
		{
			strong_enforce( boundary_conditions(ux == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(-0.3), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(0.0), {bottom_boundar_tag}) );
			
			if(is_multibody) {
				strong_enforce( boundary_conditions(ux == coeff(0.),  {3, 4, 5}) );
				strong_enforce( boundary_conditions(uy == coeff(0.),  {3, 4, 5}) );

				strong_enforce( boundary_conditions(uz == coeff(-0.08), {4}) );
				strong_enforce( boundary_conditions(uz == coeff(0.08), {5}) );
			}
		}

		std::string mesh_file;
		int top_boundary_tag;
		int bottom_boundar_tag;
		std::vector<std::pair<int, int> > contact_flags;
		bool is_multibody;
		double search_radius;
	};


	class ExampleProblem2D : public ContactProblem::ElasticityBoundaryConditions {
	public:
		ExampleProblem2D() {
			mesh_file = "../data/fine_contact_2d.e";
			top_boundary_tag = 4;
			bottom_boundar_tag = 2;
			contact_flags = {{102, 101}};
			search_radius = 0.1;
		}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy) override
		{
			strong_enforce( boundary_conditions(ux == coeff(0.),    {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(-0.15), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),    {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.0),   {bottom_boundar_tag}) );
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &, LibMeshFEFunction &) override {}

		std::string mesh_file;
		int top_boundary_tag;
		int bottom_boundar_tag;
		std::vector<std::pair<int, int> > contact_flags;
		double search_radius;
	};

	void run_contact_test(LibMeshInit &init)
	{
		auto mesh = make_shared<Mesh>(init.comm());
		ContactProblem p;
			
		//---------------------------------------------------
		// auto e_problem = make_shared<ExampleProblem2D>();
		auto e_problem = make_shared<ExampleProblem3D>();
		//---------------------------------------------------
		
		mesh->read(e_problem->mesh_file);
		p.init(init, mesh, e_problem, e_problem->contact_flags, e_problem->search_radius);
		p.step();
		p.save();
	}
}

