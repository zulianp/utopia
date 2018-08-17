#include "utopia_WearSimulation.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_LinearElasticity.hpp"
#include "utopia_NeoHookean.hpp"
#include "utopia_SaintVenantKirchoff.hpp"
#include "utopia_ContactSolver.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"

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

                    // is.read("time", [this](InputStream &is) {
                    //     is.read("dt", dt);
                    //     is.read("steps", n_time_teps);
                    // });
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
                    // double value = 0.;
                    std::string value;
                    std::string type = "volume";

                    is.read("block", block);
                    is.read("coord", coord);
                    is.read("value", value);
                    is.read("type", type);

                    if(type == "surface") {
                        auto v = test(V[coord]);
                        auto l_form = surface_integral(inner(symbolic(value), v), block);

                        auto ff = std::make_shared<ConstantForcingFunction<DVectord>>();
                        ff->init(l_form);
                        forcing_function->add(ff);
                    } else {
                        auto v = test(V[coord]);
                        auto l_form = integral(inner(symbolic(value), v), block);

                        auto ff = std::make_shared<ConstantForcingFunction<DVectord>>();
                        ff->init(l_form);
                        forcing_function->add(ff);
                    }

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

        inline Path output_path() const
        {
            return desc_.output_path;
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

                std::string type;
                is.read("type", type);

                is_steady = false;
                n_transient_steps = 1;

                if(type == "steady") {
                    is_steady = true;
                }

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
        bool is_steady;
        int n_transient_steps;

    };

    class WearSimulation::Input : public ContactSimulation {
    public:
    	Input()
        : wear_coefficient(7e-3), extrapolation_factor(10.)
    	{}

    	virtual bool init(libMesh::Parallel::Communicator &comm, InputStream &is) override
    	{
    		if(!ContactSimulation::init(comm, is)) return false;

            is.read("wear", [this](InputStream &is) {
                is.read("n-cycles",   n_cycles);
                is.read("gait-cycle", gc);
                is.read("coeff", wear_coefficient);
                is.read("extrapolation", extrapolation_factor);
            });

            gc.init(mesh->mesh_dimension());

    		return true;
    	}


    	GaitCycle gc;
        int n_cycles;
        double wear_coefficient, extrapolation_factor;


        virtual void describe(std::ostream &os) const override
        {
            ContactSimulation::describe(os);
            contact_params.describe(os);

            std::cout << "n_cycles: " << n_cycles << std::endl;
            gc.describe(os);
        }
    };

    void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &conf_file_path)
    {
        typedef utopia::ContactSolver<DSMatrixd, DVectord> ContactSolverT;
        typedef utopia::ContactStabilizedNewmark<DSMatrixd, DVectord> TransientContactSolverT;

    	Input in;
    	auto is_ptr = open_istream(conf_file_path);
    	if(!is_ptr) {
    		std::cerr << "[Error] invalid path " << conf_file_path << std::endl;
    		assert(false);
    		return;
    	}

    	in.init_sim(init.comm(), *is_ptr);
    	in.describe(std::cout);

        std::shared_ptr<ContactSolverT> solver;
        std::shared_ptr<TransientContactSolverT> transient_solver;
        if(in.is_steady) {
            solver = std::make_shared<ContactSolverT>(make_ref(in.V), in.material, in.contact_params);
        } else {
            double dt = in.gc.dt / in.n_transient_steps;
            transient_solver = std::make_shared<TransientContactSolverT>(make_ref(in.V), in.material, dt, in.contact_params);
            solver = transient_solver;
        }

        solver->set_tol(5e-6);
        solver->set_max_outer_loops(40);
        //HERE
        solver->set_use_ssn(true);

        solver->tao().atol(1e-8);
        solver->tao().rtol(1e-8);
        solver->tao().stol(1e-8);
        solver->tao().verbose(true);

        // auto ls = std::make_shared<Factorization<DSMatrixd, DVectord>>();
        // auto ls = std::make_shared<GMRES<DSMatrixd, DVectord>>();
        // ls->atol(1e-15);
        // ls->rtol(1e-15);
        // ls->stol(1e-15);
        // ls->max_it(1000);
        // // ls->verbose(true);
        // solver->set_linear_solver(ls);
        // solver->set_bypass_contact(true);
        // solver->set_max_outer_loops(10);


        Wear wear;
        wear.set_params(in.wear_coefficient, in.extrapolation_factor);
        wear.init_aux_system(
            *in.equation_systems
        );

        DVectord overriden_displacement = local_zeros(in.V.subspace(0).dof_map().n_local_dofs());
        DVectord wear_displacement = overriden_displacement;
        MechanicsState state;


        libMesh::Nemesis_IO(*in.mesh).write_timestep(in.output_path() / "wear_in.e", *in.equation_systems, (1), in.gc.t);

        libMesh::Nemesis_IO wear_io(*in.mesh);
        libMesh::Nemesis_IO wear_ovv_io(*in.mesh);
        libMesh::Nemesis_IO io(*in.mesh);

        wear_io.write_timestep(in.output_path() / "wear.e", *in.equation_systems, 1, 0);

        for(int i = 1; i <= in.n_cycles; ++i) {

            std::cout << "cycle " << i << std::endl;



            //gait-cycle
            for(int t = 0; t < in.gc.n_time_steps; ++t) {
                std::cout << "\t step " << t << std::endl;
                in.material->clear();
                in.gc.set_time_step(t);

                //set-up experiment
                    //transform mesh

                overriden_displacement.set(0.);
                in.gc.override_displacement(*in.mesh, in.V.subspace(0).dof_map(), overriden_displacement);

                {
                    auto &sys = in.equation_systems->get_system<libMesh::LinearImplicitSystem>("wear");
                    DVectord temp = overriden_displacement + wear_displacement;
                    convert(temp, *sys.solution);
                    sys.solution->close();
                    wear_ovv_io.write_timestep(in.output_path() / "wear_overr.e", *in.equation_systems, ((i-1) * in.gc.n_time_steps + t+1),((i-1) * in.gc.n_time_steps) * in.gc.dt + in.gc.t);
                }

                apply_displacement(overriden_displacement + wear_displacement, in.V.subspace(0).dof_map(), *in.mesh);

                //solve
                if(in.is_steady) {
                    solver->solve_steady();
                } else {
                    transient_solver->initial_condition(1.);
                    transient_solver->solve_dynamic(in.n_transient_steps);
                }

                //update-wear
                    //compute velocity
                    //compute normal stress
                    //compute sliding distance
                    //compute wear

                state.displacement           = solver->displacement();
                state.displacement_increment = solver->displacement();

                state.velocity               = (1./in.gc.dt) * solver->displacement();

                in.forcing_function->eval(state.displacement, state.external_force);

                in.material->stress(state.displacement, state.stress);

                const double mag_stress   = norm2(state.stress);
                const double mag_velocity = norm2(state.velocity);

                // std::cout << "mag_stress:   " << mag_stress   << std::endl;
                // std::cout << "mag_velocity: " << mag_velocity << std::endl;

                state.t = in.gc.t;

                wear.update_aux_system(
                    0,
                    state,
                    solver->contact(),
                    in.gc.dt,
                    *in.equation_systems
                );

                //clean-up experiment
                    //transform-back mesh
                apply_displacement(-(overriden_displacement + wear_displacement), in.V.subspace(0).dof_map(), *in.mesh);

                auto &sys = in.equation_systems->get_system<libMesh::LinearImplicitSystem>("wear");
                DVectord temp = state.displacement + overriden_displacement + wear_displacement;
                convert(temp, *sys.solution);
                sys.solution->close();

                io.write_timestep(in.output_path() / "gait_cycles.e", *in.equation_systems, ((i-1) * in.gc.n_time_steps + t+1),((i-1) * in.gc.n_time_steps) * in.gc.dt + in.gc.t);
            }

            //deform geometry
            wear.mesh_displacement(in.V, in.contact_surfaces, wear_displacement);
            wear_io.write_timestep(in.output_path() / "wear.e", *in.equation_systems, (i+1), i);
        }

        std::ofstream os(in.output_path() / "wear_profile.txt");

        if(os.good()) {
            wear.print(os);
        }

        os.close();

        std::cout << "finished" << std::endl;
    }

    WearSimulation::WearSimulation()
    {}

    WearSimulation::~WearSimulation()
    {}
}
