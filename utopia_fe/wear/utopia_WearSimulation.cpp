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

            equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
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


            //This should be moved
            V[0].initialize();

            forcing_function = std::make_shared<CompositeForcingFunction<DVectord>>();
            bool has_force = false;
            is.read("forcing-functions", [this, &has_force](InputStream &is) {
                is.read_all([this, &has_force](InputStream &is) {

                    int block = -1;
                    int coord = 0;
                    double value = 0.;

                    is.read("block", block);
                    is.read("coord", coord);
                    is.read("value", value);

                    auto v = test(V[coord]);
                    auto l_form = integral(inner(coeff(value), v), block);

                    auto ff = std::make_shared<ConstantForcingFunction<DVectord>>();
                    ff->init(l_form);
                    forcing_function->add(ff);

                    has_force = true;

                });
            });

            if(has_force) {
                material = std::make_shared<ForcedMaterial<DSMatrixd, DVectord>>(
                    material,
                    forcing_function
                );
            }

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
        std::shared_ptr<CompositeForcingFunction<DVectord> > forcing_function;

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


            std::set<int> temp;
    		is.read("contact", [this,&temp](InputStream &is) {
    			is.read("radius", contact_params.search_radius);

                is.read("pairs", [this,&temp](InputStream &is) {
                    is.read_all([this,&temp](InputStream &is) {
                        int master = -1, slave = -1;
                        is.read("master", master);
                        is.read("slave", slave);

                        // std::cout << master << " " << slave << std::endl;

                        assert(master != -1);
                        assert(slave  != -1);
                        temp.insert(master);
                        temp.insert(slave);

                        contact_params.contact_pair_tags.push_back({ master, slave });
                    });
                });
    		});

            contact_surfaces.clear();
            contact_surfaces.insert(contact_surfaces.end(), temp.begin(), temp.end());

            return true;
        }

        virtual void describe(std::ostream &os) const override
        {
        	ElasticitySimulation::describe(os);
        	contact_params.describe(os);
        }

        ContactParams contact_params;
        std::vector<int> contact_surfaces;

    };

    class WearSimulation::Input : public ContactSimulation {
    public:
    	Input()
    	{}

    	virtual bool init(libMesh::Parallel::Communicator &comm, InputStream &is) override
    	{
    		bool ok = ContactSimulation::init(comm, is);
            is.read("n-cycles",   n_cycles);
            is.read("gait-cycle", gc);
    		return ok;
    	}


    	GaitCycle gc;
        int n_cycles;
    };

    void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &conf_file_path)
    {
        typedef utopia::ContactSolver<DSMatrixd, DVectord> ContactSolverT;

    	Input in;
    	auto is_ptr = open_istream(conf_file_path);
    	if(!is_ptr) {
    		std::cerr << "[Error] invalid path " << conf_file_path << std::endl;
    		assert(false);
    		return;
    	}

    	in.init_sim(init.comm(), *is_ptr);
    	in.describe(std::cout);

        ContactSolverT solver(make_ref(in.V), in.material, in.contact_params);
        solver.set_tol(5e-6);

        // auto ls = std::make_shared<Factorization<DSMatrixd, DVectord>>();
        // auto ls = std::make_shared<GMRES<DSMatrixd, DVectord>>();
        // ls->atol(1e-15);
        // ls->rtol(1e-15);
        // ls->stol(1e-15);
        // ls->max_it(1000);
        // // ls->verbose(true);
        // solver.set_linear_solver(ls);
        // solver.set_bypass_contact(true);
        solver.set_max_outer_loops(10);

        libMesh::Nemesis_IO io(*in.mesh);

        Wear wear;
        wear.init_aux_system(
            *in.equation_systems
        );


        DVectord overriden_displacement = local_zeros(in.V.subspace(0).dof_map().n_local_dofs());
        MechanicsState state;

        for(int i = 0; i < in.n_cycles; ++i) {

            //gait-cycle
            for(int t = 0; t < in.gc.n_time_steps; ++t) {
                in.material->clear();

                //set-up experiment
                    //transform mesh

                in.gc.override_displacement(*in.mesh, in.V.subspace(0).dof_map(), overriden_displacement);
                apply_displacement(overriden_displacement, in.V.subspace(0).dof_map(), *in.mesh);

                //solve
                solver.solve_steady();

                //update-wear
                    //compute velocity
                    //compute normal stress
                    //compute sliding distance
                    //compute wear

                state.displacement           = solver.displacement();
                state.displacement_increment = solver.displacement();

                state.velocity               = (1./in.gc.dt) * solver.displacement();

                //TODO
                // state.external_force;//

                in.material->stress(state.displacement, state.stress);

                // state.stress -= state.external_force;
                state.t = in.gc.t;

                wear.update_aux_system(
                    0,
                    state,
                    solver.contact(),
                    in.gc.dt,
                    *in.equation_systems
                );

                //clean-up experiment
                    //transform-back mesh
                apply_displacement(-overriden_displacement, in.V.subspace(0).dof_map(), *in.mesh);

            }

            //deform geometry
        }
    }

    WearSimulation::WearSimulation()
    {}

    WearSimulation::~WearSimulation()
    {}
}
