#include "utopia_ContactTest.hpp"

#include "utopia.hpp"

//fe extension
#include "utopia_fe_core.hpp"
#include "MortarAssembler.hpp"
#include "ParMortarAssembler.hpp"
#include "utopia_Socket.hpp"
#include "utopia_ContactSimParams.hpp"
#include "utopia_Polygon.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"
#include "utopia_ContactProblem.hpp"
#include "utopia_assemble_contact.hpp"

#include "moonolith_profiler.hpp"

#include <libmesh/mesh.h>
#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/parallel_mesh.h>
#include "libmesh/linear_partitioner.h"
#include "LibmeshContactForMoose.hpp"
#include "LibmeshTransferForMoose.hpp"
#include "LibmeshTransferForMooseReverse.hpp"
#include <libmesh/mesh_base.h>


#include <iostream>

using namespace libMesh;
using std::make_shared;
using std::shared_ptr;

namespace utopia {

	class ExampleProblemBase : public ContactProblem::ElasticityBoundaryConditions {
	public:
		virtual ~ExampleProblemBase() {}
		std::string mesh_file;
		int top_boundary_tag;
		int bottom_boundar_tag;
		std::vector<std::pair<int, int> > contact_flags;
		double search_radius;
		double dt;
		int n_steps;
	};


	class ExampleProblem3D : public ExampleProblemBase {
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
			strong_enforce( boundary_conditions(uz == coeff(-dt*0.3), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(0.0), {bottom_boundar_tag}) );
			
			if(is_multibody) {
				strong_enforce( boundary_conditions(ux == coeff(0.),  {3, 4, 5}) );
				strong_enforce( boundary_conditions(uy == coeff(0.),  {3, 4, 5}) );

				strong_enforce( boundary_conditions(uz == coeff(-dt*0.08), {4}) );
				strong_enforce( boundary_conditions(uz == coeff(dt*0.08), {5}) );
			}
		}

		bool is_multibody;
	};

	class QuasiHertz : public ExampleProblemBase {
	public:
		QuasiHertz() {
			mesh_file = "../data/quasi_hertz_8732.e";
			top_boundary_tag = 2;
			bottom_boundar_tag = 4;
			contact_flags = {{1, 3}};
			search_radius = 0.4;
			double dt = 1.;
			n_steps = 1;
		}

		void set_up_fine()
		{
			mesh_file = "../data/quasi_hertz_260310.e";
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &)  override {}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) override
		{
			strong_enforce( boundary_conditions(ux == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(-dt*0.1), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(0.0), {bottom_boundar_tag}) );
		}
	};

	class QuasiSignorini : public ExampleProblemBase {
	public:
		QuasiSignorini() {
			mesh_file = "../data/quasi_signorini_4593.e";
			top_boundary_tag = 2;
			bottom_boundar_tag = 4;
			contact_flags = {{3, 1}};
			search_radius = 0.2;
			double dt = 1.;
			n_steps = 1;
		}

		void set_up_middle_res()
		{
			mesh_file = "../data/quasi_signorini_32894.e";
		}

		void set_up_adaptive()
		{
			mesh_file = "../data/half_sphere_with_plate.e";
		}

		void set_up_fine_res()
		{
			mesh_file = "../data/quasi_signorini_1928354.e";
		}

		void set_up_time_dependent(){
			dt = 0.1;
			search_radius *= dt;
			n_steps = 30;
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &)  override {}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) override
		{
			strong_enforce( boundary_conditions(ux == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(-dt*0.1), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(0.0), {bottom_boundar_tag}) );
		}
	};


	class ExampleProblem2D : public ExampleProblemBase {
	public:
		ExampleProblem2D() {
			mesh_file = "../data/fine_contact_2d.e";
			top_boundary_tag = 4;
			bottom_boundar_tag = 2;
			contact_flags = {{102, 101}};
			search_radius = 0.1;
			dt = 1.;
			n_steps = 1;
		}

		void set_up_coarse_t()
		{
			mesh_file = "../data/contact_2d.e";
			search_radius = 0.05;
			contact_flags = {{101, 102}};
			dt = 0.01;
			n_steps = 100;
		}

		void set_up_m_coarse_t()
		{
			// mesh_file = "../data/m_contact_2d.e";
			mesh_file = "../data/m_contact_2d_ref.e";
			search_radius = 0.01;
			// contact_flags = {{101, 102}, {103, 104}, {105, 106}};
			contact_flags = {{102, 101}, {103, 104}, {105, 106}};
			dt = 0.025;
			n_steps = 100;
		}


		void set_up_coarse()
		{
			mesh_file = "../data/contact_2d.e";
			search_radius = 0.5;
			contact_flags = {{101, 102}};
			dt = 1.;
			n_steps = 1;
		}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy) override
		{
			strong_enforce( boundary_conditions(ux == coeff(0.),    	 {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(-dt * 0.15), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),    	 {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.0),   	 {bottom_boundar_tag}) );
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &, LibMeshFEFunction &) override {}
	};

	void run_contact_test(LibMeshInit &init)
	{
		auto mesh = make_shared<Mesh>(init.comm());
		ContactProblem p;
			
		//---------------------------------------------------
		// auto e_problem = make_shared<ExampleProblem2D>();
		// e_problem->set_up_m_coarse_t();
		
		// auto e_problem = make_shared<ExampleProblem3D>();
		// auto e_problem = make_shared<QuasiHertz>();
		auto e_problem = make_shared<QuasiSignorini>(); 
		// e_problem->set_up_fine_res();
		// e_problem->set_up_adaptive();
		// e_problem->set_up_time_dependent();
		//---------------------------------------------------
		
		std::cout << "reading mesh...." << std::flush;
		// mesh->set_distributed();
		mesh->read(e_problem->mesh_file);
		// std::cout << "is_replicated: " << mesh->is_replicated() << std::endl;
		mesh->delete_remote_elements();
		
		std::cout << "done" << std::endl;

		p.init(init, mesh, e_problem, e_problem->contact_flags, e_problem->search_radius);
		

		double t = 0.0;
		for(int i = 0; i < e_problem->n_steps; ++i) {
			std::cout << "t: " << t << "/" << (e_problem->dt * e_problem->n_steps) << std::endl;
			std::cout << std::flush;
			t += e_problem->dt;
			p.step();			
		}

		p.save();
	}
}

