#include "utopia_WearSimulation.hpp"

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml_print.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_LinearElasticity.hpp"
#include "utopia_NeoHookean.hpp"
#include "utopia_SaintVenantKirchoff.hpp"
#include "utopia_ContactSolver.hpp"

#include <iostream>
#include <fstream>

namespace utopia {

	template class ContactSolver<DSMatrixd, DVectord>;
	typedef utopia::ContactSolver<DSMatrixd, DVectord> ContactSolverT;
	typedef std::array<double, 2> Point2d;
	typedef std::array<double, 3> Point3d;
	typedef std::function<std::array<double, 2>(const std::array<double, 2> &p)> Fun2d;
	typedef std::function<std::array<double, 3>(const std::array<double, 3> &p)> Fun3d;


	class WearSimulation::GaitCycle {
	public:
		void init(
			const int dim,
			rapidxml::xml_node<> &xml_gait_cycle)
		{
			using namespace rapidxml;
		}

		double start_angle_radian;
		double angle_degree;
		double angle_radian;
		double d_angle;
		Fun2d rotate2;
		Fun3d rotate3;
	};

	class WearSimulation::Input {
	public:

		Input()
		: output_path("./")
		{}

		bool init_from_xml(libMesh::LibMeshInit &init, std::istream &is)
		{
			using namespace rapidxml;

			file<> f(is);
			xml_document<> doc;  
			doc.parse<0>(f.data());  

			xml_node<> *wear_simulation = doc.first_node("simulation");
			
			if(!wear_simulation) {
				std::cerr << "could not find <simulation> node" << std::endl;
				return false;
			}

			xml_node<> *xml_mesh = wear_simulation->first_node("mesh");
			if(!xml_mesh) {
				std::cerr << "could not find <mesh> node" << std::endl;
				return false;
			}

			std::cout << "mesh:   " << xml_mesh->value() << std::endl;

			mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
			mesh->read(xml_mesh->value());
			auto dim = mesh->mesh_dimension();

			auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
			auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("wear");
			main_sys_num = sys.number();

			auto elem_order = libMesh::FIRST;

			auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
			auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
			auto V  = Vx * Vy;

			if(dim == 3) {
				V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
			}

			/////////////////////////////////////////////////////////////////////////////

			xml_node<> *xml_bc = wear_simulation->first_node("boundary-conditions");
			xml_node<> *xml_dirichlet = xml_bc->first_node("dirichlet");

			std::cout << "side\tcoord\tvalue" << std::endl;
			for (xml_node<> *xml_cond = xml_dirichlet->first_node(); xml_cond; xml_cond = xml_cond->next_sibling())
			{
				xml_attribute<> *cond_side_attr  = xml_cond->first_attribute("side"); assert(cond_side_attr);
				xml_attribute<> *cond_coord_attr = xml_cond->first_attribute("coord"); assert(cond_coord_attr);

				auto side_set = atoi(cond_side_attr->value());
				auto coord    = atoi(cond_coord_attr->value());
				auto value    = atof(xml_cond->value());

				std::cout << side_set << "\t" << coord << "\t" << value << std::endl;

				auto u = trial(V[coord]);
				init_constraints(constraints(
					boundary_conditions(u == coeff(value), {side_set})
				));
			}

			/////////////////////////////////////////////////////////////////////////////

			xml_node<> *xml_model 	   = wear_simulation->first_node("model");
			xml_node<> *xml_material   = xml_model->first_node("material");
			xml_node<> *xml_parameters = xml_model->first_node("parameters");
			
			if(xml_parameters) {
				xml_node<> *xml_mu     = xml_parameters->first_node("mu");
				xml_node<> *xml_lambda = xml_parameters->first_node("lambda");

				const double mu     = xml_mu?     atof(xml_mu->value()) : 1.;
				const double lambda = xml_lambda? atof(xml_lambda->value()) : 1.;

				params = LameeParameters(mu, lambda);
			}

			if(xml_material) {
				if(xml_material->value() == std::string("NeoHookean")) {
					std::cout << "material: NeoHookean" << std::endl;
					material = std::make_shared<NeoHookean<decltype(V), DSMatrixd, DVectord>>(V, params);
				} else if(xml_model->value() == std::string("SaintVenantKirchoff")) {
					std::cout << "material: SaintVenantKirchoff" << std::endl;
					material = std::make_shared<SaintVenantKirchoff<decltype(V), DSMatrixd, DVectord>>(V, params);
				} else //if(mode->value() == "LinearElasticity") 
				{
					std::cout << "material: LinearElasticity" << std::endl;
					material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, params);
				}
			} else {
				material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, params);
			}

			params.describe(std::cout);

			/////////////////////////////////////////////////////////////////////////////////////
			xml_node<> *xml_contact = wear_simulation->first_node("contact"); assert(xml_contact);
			xml_attribute<> *radius_att = xml_contact->first_attribute("radius");

			if(radius_att) {
				contact_params.search_radius = atof(radius_att->value());
			}

			std::cout << "master\tslave" << std::endl;
			for (xml_node<> *xml_pair = xml_contact->first_node(); xml_pair; xml_pair = xml_pair->next_sibling())
			{
				xml_node<> *xml_master = xml_pair->first_node("master"); assert(xml_master);
				xml_node<> *xml_slave  = xml_pair->first_node("slave");  assert(xml_slave);
				
				contact_params.contact_pair_tags.push_back({
					atoi(xml_master->value()),
					atoi(xml_slave->value()),
				});

				std::cout << contact_params.contact_pair_tags.back().first << "\t"	
						  << contact_params.contact_pair_tags.back().second << std::endl;
			}

			/////////////////////////////////////////////////////////////////////////////////////

			xml_node<> *xml_time  = wear_simulation->first_node("time");
			xml_node<> *xml_dt    = xml_time->first_node("dt");
			xml_node<> *xml_steps = xml_time->first_node("steps"); 

			dt = atof(xml_dt->value());
			n_time_teps = atoi(xml_steps->value());

			std::cout << "dt: " << dt << " n_time_teps: " << n_time_teps << std::endl;

			/////////////////////////////////////////////////////////////////////////////////////

			auto xml_output = wear_simulation->first_node("ouput");
			if(xml_output) {
				output_path = xml_output->value();
			}

			/////////////////////////////////////////////////////////////////////////////////////
			return true;
		}

		std::string mesh_path;
		std::string output_path;
		LameeParameters params;
		std::shared_ptr<libMesh::DistributedMesh> mesh;
		ProductFunctionSpace<LibMeshFunctionSpace> V;
		std::shared_ptr<libMesh::EquationSystems> equation_systems;
		std::shared_ptr<ElasticMaterial<DSMatrixd, DVectord> > material;
		ContactParams contact_params;
		int main_sys_num;
		int aux_sys_num;
		int n_time_teps;
		double dt;
	};

	void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &xml_file_path)
	{
		Input in;
		std::ifstream is(xml_file_path);
		in.init_from_xml(init, is);
	}

	WearSimulation::WearSimulation()
	{}

	WearSimulation::~WearSimulation()
	{}
}
