#include "utopia_WearSimulation.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_LinearElasticity.hpp"
#include "utopia_NeoHookean.hpp"
#include "utopia_SaintVenantKirchoff.hpp"
#include "utopia_ContactSolver.hpp"

#include "utopia_ui.hpp"
#include "utopia_GaitCycle.hpp"
#include "utopia_Wear.hpp"

#include <iostream>
#include <fstream>

namespace utopia {

	class ElasticitySimulation {
	public:
		virtual ~ElasticitySimulation() {}

		ElasticitySimulation() :
		mesh_path("../data/mesh2.e"),
		output_path("output2"),
		material_name("LinearElasticity"),
		params(1, 1),
		n_time_teps(10),
		dt(0.1)
		{}

		virtual bool init_sim(libMesh::Parallel::Communicator &comm, InputStream &is)
		{
			bool ok = false;
			if(is.object_begin("simulation"))
			{
				ok = init(comm, is);
			}

			is.object_end();
			return ok;
		}

		virtual bool init(libMesh::Parallel::Communicator &comm, InputStream &is)
		{
			is.read("mesh", mesh_path);

			mesh = std::make_shared<libMesh::DistributedMesh>(comm);
			mesh->read(mesh_path);
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

			if(is.object_begin("model")) {
				is.read("material", material_name);

				std::cout << material_name << std::endl;

				if(material_name == "NeoHookean") {
					material = std::make_shared<NeoHookean<decltype(V), DSMatrixd, DVectord>>(V, params);
				} else if(material_name == "SaintVenantKirchoff") {
					material = std::make_shared<SaintVenantKirchoff<decltype(V), DSMatrixd, DVectord>>(V, params);
                } else /*if(material_name == "LinearElasticity")*/ {
					material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, params);
				}

				if(is.object_begin("parameters")) {
					double mu, lambda;
					is.read("mu", mu);
					is.read("lambda", lambda);
					params = LameeParameters(mu, lambda);
				}

                is.object_end(); //parameters
            }

            is.object_end(); //model

            /////////////////////////////////////////////////////////////////////////////

            if(is.object_begin("boundary-conditions")) {

            	if(is.object_begin("dirichlet")) {

            		is.start();

            		while(is.good()) {
            			int side_set = 0, coord = 0;
            			double value = 0;

            			is.read("side", side_set);
            			is.read("coord", coord);
            			is.read("value", value);

            			std::cout << side_set << ", " << coord << ", " << value << std::endl;

            			auto u = trial(V[coord]);
            			init_constraints(constraints(
            				boundary_conditions(u == coeff(value), {side_set})
            			));

            			is.next();
            		}

            		is.finish();
            	}

                is.object_end(); //dirichlet
            }

            is.object_end(); //boundary-conditions

            /////////////////////////////////////////////////////////////////////////////

            if(is.object_begin("time")) {
            	is.read("dt", dt);
            	is.read("steps", n_time_teps);
            }

            is.object_end(); //time

            /////////////////////////////////////////////////////////////////////////////

            is.read("output", output_path);


            if(!material) {
            	material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, params);
            }

            return true;
        }

        virtual void describe(std::ostream &os) const
        {
        	os << "mesh_path:\t" << mesh_path << "\n";
        	os << "output_path:\t" << output_path << "\n";
        	os << "material_name:\t" << material_name << "\n";
        	os << "material_params:\n";
        	params.describe(os);
        	os << "n_time_teps:\t" << n_time_teps << "\n";
        	os << "dt:\t" << dt << "\n";
        }

    public:
    	std::shared_ptr<libMesh::DistributedMesh> mesh;
    	ProductFunctionSpace<LibMeshFunctionSpace> V;
    	std::shared_ptr<libMesh::EquationSystems> equation_systems;
    	std::shared_ptr<ElasticMaterial<DSMatrixd, DVectord> > material;

    	int main_sys_num;
    	int aux_sys_num;

    private:
    	std::string mesh_path;
    	std::string output_path;
    	std::string material_name;
    	LameeParameters params;
    	int n_time_teps;
    	double dt;
    };

    class ContactSimulation : public ElasticitySimulation {
    public:
    	virtual ~ContactSimulation() {}

    	virtual bool init(libMesh::Parallel::Communicator &comm, InputStream &is) override
    	{
    		bool ok = true;
    		if(!ElasticitySimulation::init(comm, is)) {
    			ok = false;
    		}

    		if(is.object_begin("contact")) {
    			is.read("radius", contact_params.search_radius);

    			is.start("pair");

    			while(is.good()) {
    				int master = -1, slave = -1;
    				is.read("master", master);
    				is.read("slave", slave);

    				assert(master != -1);
    				assert(slave != -1);

    				contact_params.contact_pair_tags.push_back({ master, slave });

    				is.next();
    			}

    			is.finish();
    		}

            is.object_end(); //contact
            return true;
        }

        virtual void describe(std::ostream &os) const override
        {
        	ElasticitySimulation::describe(os);
        	contact_params.describe(os);
        }

        ContactParams contact_params;

    };

    class WearSimulation::Input : public ContactSimulation {
    public:
    	Input()
    	{}

    	virtual bool init(libMesh::Parallel::Communicator &comm, InputStream &is) override
    	{
    		bool ok = ContactSimulation::init(comm, is);

    		if(is.object_begin("gait-cycle")) {
    			gc.init(mesh->mesh_dimension(), is);
    		} else {
    			ok = false;
    		}

    		is.object_end();
    		return ok;
    	}

    	GaitCycle gc;
    };

    void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &conf_file_path)
    {
    	Input in;
    	auto is_ptr = open_istream(conf_file_path);
    	if(!is_ptr) {
    		std::cerr << "[Error] invalid path " << conf_file_path << std::endl;
    		assert(false);
    		return;

    	}

    	in.init_sim(init.comm(), *is_ptr);
    	in.describe(std::cout);
    }

    WearSimulation::WearSimulation()
    {}

    WearSimulation::~WearSimulation()
    {}
}
