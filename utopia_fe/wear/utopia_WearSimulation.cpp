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

        class Desc : public Serializable {
        public:
            Desc() :
                mesh_path("../data/mesh2.e"),
                output_path("output2"),
                material_name("LinearElasticity"),
                params(1, 1),
                n_time_teps(10),
                dt(0.1)
            { }

            void read(InputStream &is) {
               is.read("mesh", mesh_path);

               is.read("model", [this](InputStream &is) {
                    is.read("material", material_name);

                    is.read("parameters", [this](InputStream &is) {
                        is.read("mu", params.default_mu);
                        is.read("lambda", params.default_lambda);
                    });

                    is.read("time", [this](InputStream &is) {
                        is.read("dt", dt);
                        is.read("steps", n_time_teps);
                    });
               });

               is.read("output", output_path);
            }

            std::string mesh_path;
            std::string output_path;
            std::string material_name;
            LameeParameters params;
            int n_time_teps;
            double dt;
        };

		ElasticitySimulation()
		{}

		virtual bool init_sim(libMesh::Parallel::Communicator &comm, InputStream &is)
		{
			bool ok = false;

			is.read("simulation", [this, &ok, &comm](InputStream &is) {
                ok = init(comm, is);
            });

			return ok;
		}

		virtual bool init(libMesh::Parallel::Communicator &comm, InputStream &is)
		{
            is.read(desc_);

			mesh = std::make_shared<libMesh::DistributedMesh>(comm);
			mesh->read(desc_.mesh_path);
			auto dim = mesh->mesh_dimension();

			auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
			auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("wear");
			main_sys_num = sys.number();

			auto elem_order = libMesh::FIRST;

			auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
			auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
			V  = Vx * Vy;

			if(dim == 3) {
				V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
			}

            /////////////////////////////////////////////////////////////////////////////

			if(desc_.material_name == "NeoHookean") {
				material = std::make_shared<NeoHookean<decltype(V), DSMatrixd, DVectord>>(V, desc_.params);
			} else if(desc_.material_name == "SaintVenantKirchoff") {
				material = std::make_shared<SaintVenantKirchoff<decltype(V), DSMatrixd, DVectord>>(V, desc_.params);
            } else /*if(desc_.material_name == "LinearElasticity")*/ {
				material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, desc_.params);
			}

            /////////////////////////////////////////////////////////////////////////////

            is.read("boundary-conditions", [this](InputStream &is) {
            	is.read("dirichlet", [this](InputStream &is) {

                    is.read_all([this](InputStream &is) {
                        int side_set = 0, coord = 0;
                        double value = 0;
                        std::string expr = "0";

                        is.read("side", side_set);
                        is.read("coord", coord);
                        is.read("value", expr);

                        // std::cout << side_set << " " << coord << " " << value << std::endl;

                        auto u = trial(V[coord]);
                        init_constraints(constraints(
                            // boundary_conditions(u == coeff(value), {side_set})
                            boundary_conditions(u == symbolic(expr), {side_set})
                        ));

                    });

            	});
            });

            /////////////////////////////////////////////////////////////////////////////
            return true;
        }

        virtual void describe(std::ostream &os) const
        {
        	os << "mesh_path:\t" << desc_.mesh_path << "\n";
        	os << "output_path:\t" << desc_.output_path << "\n";
        	os << "material_name:\t" << desc_.material_name << "\n";
        	os << "material_params:\n";
        	desc_.params.describe(os);
        	os << "n_time_teps:\t" << desc_.n_time_teps << "\n";
        	os << "dt:\t" << desc_.dt << "\n";
        }

    public:
    	std::shared_ptr<libMesh::DistributedMesh> mesh;
    	ProductFunctionSpace<LibMeshFunctionSpace> V;
    	std::shared_ptr<libMesh::EquationSystems> equation_systems;
    	std::shared_ptr<ElasticMaterial<DSMatrixd, DVectord> > material;

    	int main_sys_num;
    	int aux_sys_num;

        Desc desc_;
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

    		is.read("contact", [this](InputStream &is) {
    			is.read("radius", contact_params.search_radius);

                is.read("pairs", [this](InputStream &is) {
                    is.read_all([this](InputStream &is) {
                        int master = -1, slave = -1;
                        is.read("master", master);
                        is.read("slave", slave);

                        // std::cout << master << " " << slave << std::endl;

                        assert(master != -1);
                        assert(slave  != -1);

                        contact_params.contact_pair_tags.push_back({ master, slave });
                    });
                });
    		});

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
            is.read("gait-cycle", gc);
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
