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

#include <iostream>
#include <fstream>

namespace utopia {
	class WearSimulation::Impl {
	public:

		Impl()
		: output_path("./")
		{}

		bool run(libMesh::LibMeshInit &init, std::istream &is)
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

			for (xml_node<> *xml_cond = xml_dirichlet->first_node(); xml_cond; xml_cond = xml_cond->next_sibling())
			{
				xml_attribute<> *cond_side_attr  = xml_cond->first_attribute("side"); assert(cond_side_attr);
				xml_attribute<> *cond_coord_attr = xml_cond->first_attribute("coord"); assert(cond_coord_attr);

				auto side_set = atoi(cond_side_attr->value());
				auto coord    = atoi(cond_coord_attr->value());
				auto value    = atof(xml_cond->value());

				std::cout << side_set << " " << coord << " " << value << std::endl;

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

			/////////////////////////////////////////////////////////////////////////////////////

			
			params.describe(std::cout);
			return true;
		}

		std::string mesh_path;
		std::string output_path;
		LameeParameters params;
		std::shared_ptr<libMesh::DistributedMesh> mesh;
		ProductFunctionSpace<LibMeshFunctionSpace> V;
		std::shared_ptr<libMesh::EquationSystems> equation_systems;
		std::shared_ptr<ElasticMaterial<DSMatrixd, DVectord> > material;
	};

	void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &xml_file_path)
	{
		Impl impl;

		std::ifstream is(xml_file_path);
		impl.run(init, is);
	}

	WearSimulation::WearSimulation()
	{

	}

	WearSimulation::~WearSimulation()
	{

	}
}
