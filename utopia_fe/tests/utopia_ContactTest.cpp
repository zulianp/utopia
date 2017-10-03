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

	class ExampleProblemBase : public ContactProblem::ElasticityBoundaryConditions, public ContactProblem::ElasticityForcingFunction {
	public:
		virtual ~ExampleProblemBase() {}
		std::string mesh_file;
		int top_boundary_tag;
		int bottom_boundar_tag;
		std::vector<std::pair<int, int> > contact_flags;
		double search_radius;
		double dt;
		int n_steps;
		bool dynamic_contact;
		bool is_inpulse;

		ExampleProblemBase()
		{
			is_inpulse = false;
		}
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
			dt = 1.;
			dynamic_contact = false;

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
			dt = 1.;
			n_steps = 1;
			dynamic_contact = false;
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
			strong_enforce( boundary_conditions(uz == coeff(-0.1), {top_boundary_tag}) );

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
			dt = 1.;
			n_steps = 1;
			dynamic_contact = false;
		}

		int block_id() const override
		{
			return 1;
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
			dt = 0.3;
			search_radius *= dt;
			n_steps = 10;
		}

		void set_up_time_dependent_dynamic(){
			set_up_adaptive();
			set_up_time_dependent();
			dynamic_contact = true;
			n_steps = 100;
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &)  override {}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) override
		{
			if(!dynamic_contact) {
				strong_enforce( boundary_conditions(ux == coeff(0.), {top_boundary_tag}) );
				strong_enforce( boundary_conditions(uy == coeff(0.), {top_boundary_tag}) );
				strong_enforce( boundary_conditions(uz == coeff(-0.2), {top_boundary_tag}) );
			}

			strong_enforce( boundary_conditions(ux == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(0.0), {bottom_boundar_tag}) );
		}

		virtual void fill(libMesh::DenseVector<libMesh::Real> &v) override
		{
			if(dynamic_contact) {
				v(2) = -0.1;
			} else {
				v.zero();
			}
		}
	};

	class Rocks : public ExampleProblemBase {
	public:
		Rocks() {
			mesh_file = "../data/two_rocks_2mm_8041.e";
			top_boundary_tag = 2;
			bottom_boundar_tag = 4;
			contact_flags = {{3, 1}};
			search_radius = 2.;
			dt = 1.;
			n_steps = 9;
			dynamic_contact = false;
		}

		int block_id() const override
		{
			return 1;
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &)  override {}

		void apply(LibMeshFEFunction &ux, LibMeshFEFunction &uy, LibMeshFEFunction &uz) override
		{
			strong_enforce( boundary_conditions(ux == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.), {top_boundary_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(-0.5), {top_boundary_tag}) );

			strong_enforce( boundary_conditions(ux == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.),  {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uz == coeff(0.0), {bottom_boundar_tag}) );
		}

		virtual void fill(libMesh::DenseVector<libMesh::Real> &v) override
		{
			v.zero();
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
			dynamic_contact = false;

		}

		void set_up_dynamic()
		{
			mesh_file = "../data/coarse_contact_2d.e";
			search_radius = 0.002;
			dt = 0.1;
			dynamic_contact = true;
			n_steps = 200;
		}

		void set_up_dynamic_with_impulse()
		{
			mesh_file = "../data/coarse_contact_2d.e";
			search_radius = 0.002;
			dt = 0.02;
			dynamic_contact = true;
			n_steps = 100;
			is_inpulse = true;
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
			search_radius = 0.1;
			contact_flags = {{101, 102}, {103, 104}, {105, 106}};
			// contact_flags = {{102, 101}, {103, 104}, {105, 106}};
			dt = 0.05;
			n_steps = 60;
		}

		void set_up_m_coarse_t_dynamic()
		{
			set_up_m_coarse_t();
			dynamic_contact = true;
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
			if(!dynamic_contact) {
				strong_enforce( boundary_conditions(ux == coeff(0.),    	 {top_boundary_tag}) );
				strong_enforce( boundary_conditions(uy == coeff(-0.15), 	 {top_boundary_tag}) );
			}

			strong_enforce( boundary_conditions(ux == coeff(0.),    	 {bottom_boundar_tag}) );
			strong_enforce( boundary_conditions(uy == coeff(0.0),   	 {bottom_boundar_tag}) );
		}

		virtual void fill(libMesh::DenseVector<libMesh::Real> &v) override
		{
			if(dynamic_contact) {
				if(is_inpulse) {
					v(1) = -25;
				} else {
					v(1) = -1.;
				}
			} else {
				v.zero();
			}
		}

		void apply(LibMeshFEFunction &, LibMeshFEFunction &, LibMeshFEFunction &) override {}
	};

	void run_contact_test(LibMeshInit &init)
	{
		auto mesh = make_shared<DistributedMesh>(init.comm());
		// mesh->partitioner().reset(new LinearPartitioner());
		// auto mesh = make_shared<Mesh>(init.comm());
		ContactProblem p;
			
		// ---------------------------------------------------
		auto e_problem = make_shared<ExampleProblem2D>();
		e_problem->set_up_m_coarse_t_dynamic();

		// e_problem->set_up_m_coarse_t();
		// auto e_problem = make_shared<Rocks>();
		// e_problem->set_up_m_coarse_t_dynamic();
		// e_problem->set_up_dynamic_with_impulse();
		 
		
		// auto e_problem = make_shared<ExampleProblem3D>();
		// auto e_problem = make_shared<QuasiHertz>();
		
		// auto e_problem = make_shared<QuasiSignorini>(); 
		// e_problem->set_up_time_dependent_dynamic();


		// e_problem->set_up_fine_res();
		// e_problem->set_up_adaptive();
		// e_problem->set_up_time_dependent();
		
		//---------------------------------------------------
		
		MOONOLITH_EVENT_BEGIN("read_mesh");
		mesh->read(e_problem->mesh_file);
		MOONOLITH_EVENT_END("read_mesh");

		p.init(init, mesh, e_problem, e_problem, e_problem->contact_flags, e_problem->search_radius);
		p.set_dynamic_contact(e_problem->dynamic_contact);
		p.is_inpulse(e_problem->is_inpulse);
		p.save(e_problem->dt);

		double t = 0.0;
		for(int i = 0; i < e_problem->n_steps; ++i) {
			std::cout << "t: " << t << "/" << (e_problem->dt * e_problem->n_steps) << std::endl;
			std::cout << std::flush;
			t += e_problem->dt;
			p.step(e_problem->dt);	
			p.save(e_problem->dt);		
		}

		p.save_energy("energy.txt");
	}
}

