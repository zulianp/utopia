#include "utopia_WearSimulation.hpp"

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml_print.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_LinearElasticity.hpp"

#include <iostream>
#include <fstream>

namespace utopia {
	class WearSimulation::Impl {
	public:

		Impl()
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


			if(xml_model) {
				if(xml_model->value() == std::string("NeoHookean")) {

				} else if(xml_model->value() == std::string("SaintVenantKirchoff")) {

				} else //if(mode->value() == "LinearElasticity") 
				{
					material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, params);
				}
			}

			std::cout << "mesh:   " << xml_mesh->value() << std::endl;
			params.describe(std::cout);
			return true;
		}

		std::string mesh_path;
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
