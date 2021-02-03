// #include "utopia_libmesh.hpp"
// #include "utopia_WearSimulation.hpp"

// #include "utopia_fe_base.hpp"

// #ifndef UTOPIA_WITH_TRILINOS_ALGEBRA

// #include "utopia_LameeParameters.hpp"
// #include "utopia_ElasticMaterial.hpp"
// #include "utopia_LinearElasticity.hpp"
// #include "utopia_NeoHookean.hpp"
// #include "utopia_SaintVenantKirchoff.hpp"
// #include "utopia_ContactSolver.hpp"
// #include "utopia_ContactStabilizedNewmark.hpp"
// #include "utopia_UIForcingFunction.hpp"

// #include "utopia_ui.hpp"
// #include "utopia_GaitCycle.hpp"
// #include "utopia_Wear.hpp"

// #include "libmesh/mesh_refinement.h"

// #include <iostream>
// #include <fstream>

// namespace utopia {
//     static void refine(const int n_refs, libMesh::MeshBase &mesh)
//     {
//         if(n_refs <= 0) return;

//         libMesh::MeshRefinement mesh_refinement(mesh);
//         mesh_refinement.make_flags_parallel_consistent();
//         mesh_refinement.uniformly_refine(n_refs);
//     }

//     class ElasticitySimulation {
//     public:
//         virtual ~ElasticitySimulation() {}

//         class Desc : public Configurable {
//         public:
//             Desc() :
//                 mesh_path("../data/mesh2.e"),
//                 output_path("output2"),
//                 material_name("LinearElasticity"),
//                 params(1, 1),
//                 n_time_teps(10),
//                 dt(0.1)
//             { }

//             void read(Input &is) {
//                is.get("mesh", mesh_path);

//                mesh_refinements = 0;
//                is.get("mesh-refinements", mesh_refinements);

//                is.get("model", [this](Input &is) {
//                     is.get("material", material_name);

//                     is.get("parameters", [this](Input &is) {
//                         is.get("mu", params.default_mu);
//                         is.get("lambda", params.default_lambda);
//                     });

//                     stabilization = "none";

//                     is.get("stabilization", stabilization);
//                     stabilization_mag = 0.1;
//                     is.get("stabilization-mag", stabilization_mag);

//                     // is.get("time", [this](Input &is) {
//                     //     is.get("dt", dt);
//                     //     is.get("steps", n_time_teps);
//                     // });
//                });

//                is.get("output", output_path);
//             }

//             std::string mesh_path;
//             std::string output_path;
//             std::string material_name;
//             LameeParameters params;
//             int n_time_teps;
//             double dt;
//             std::string stabilization;
//             double stabilization_mag;
//             int mesh_refinements;
//         };

//         ElasticitySimulation()
//         {}

//         virtual bool init_sim(libMesh::Parallel::Communicator &comm, Input &is)
//         {
//             bool ok = false;

//             is.get("simulation", [this, &ok, &comm](Input &is) {
//                 ok = init(comm, is);
//             });

//             return ok;
//         }

//         virtual bool init(libMesh::Parallel::Communicator &comm, Input &is)
//         {
//             // is.read(desc_);
//             desc_.read(is);

//             mesh = std::make_shared<libMesh::DistributedMesh>(comm);
//             mesh->read(desc_.mesh_path);

//             refine(desc_.mesh_refinements, *mesh);

//             auto dim = mesh->mesh_dimension();

//             equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
//             auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("wear");
//             main_sys_num = sys.number();

//             auto elem_order = libMesh::FIRST;

//             auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
//             auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
//             V  = Vx * Vy;

//             if(dim == 3) {
//                 V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
//             }

//             /////////////////////////////////////////////////////////////////////////////

//             if(desc_.material_name == "NeoHookean") {
//                 material = std::make_shared<NeoHookean<decltype(V), USparseMatrix, UVector>>(V, desc_.params);
//             } else if(desc_.material_name == "SaintVenantKirchoff") {
//                 material = std::make_shared<SaintVenantKirchoff<decltype(V), USparseMatrix, UVector>>(V,
//                 desc_.params);
//             } else /*if(desc_.material_name == "LinearElasticity")*/ {
//                 material = std::make_shared<LinearElasticity<decltype(V), USparseMatrix, UVector>>(V, desc_.params);
//             }

//             if(desc_.stabilization != "none") {
//                 std::cout << "using stabilization: " << desc_.stabilization << " mag: " << desc_.stabilization_mag <<
//                 std::endl; material = std::make_shared<StabilizedMaterial<decltype(V), USparseMatrix, UVector>>(V,
//                 desc_.stabilization_mag, material, desc_.stabilization);
//             }

//             /////////////////////////////////////////////////////////////////////////////

//             is.get("boundary-conditions", [this](Input &is) {
//                 is.get("dirichlet", [this](Input &is) {

//                     is.get_all([this](Input &is) {
//                         int side_set = 0, coord = 0;

//                         is.get("side", side_set);
//                         is.get("coord", coord);

//                         auto u = trial(V[coord]);

// #ifdef UTOPIA_WITH_TINY_EXPR
//                         std::string expr = "0";
//                         is.get("value", expr);
//                         auto g = symbolic(expr);
// #else
//                         double value = 0;
//                         is.get("value", value);
//                         auto g = coeff(value);
// #endif //UTOPIA_WITH_TINY_EXPR

//                         init_constraints(constraints(
//                             boundary_conditions(u == g, {side_set})
//                         ));

//                     });

//                 });
//             });

//             //This should be moved
//             V[0].initialize();

//             auto ff = std::make_shared<UIForcingFunction<decltype(V), UVector>>(V);
//             is.get("forcing-functions", *ff);
//             forcing_function = ff;

//             material = std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(
//                 material,
//                 forcing_function
//             );

//             /////////////////////////////////////////////////////////////////////////////
//             return true;
//         }

//         virtual void describe(std::ostream &os) const
//         {
//             os << "mesh_path:\t" << desc_.mesh_path << "\n";
//             os << "output_path:\t" << desc_.output_path << "\n";
//             os << "material_name:\t" << desc_.material_name << "\n";
//             os << "material_params:\n";
//             desc_.params.describe(os);
//             os << "n_time_teps:\t" << desc_.n_time_teps << "\n";
//             os << "dt:\t" << desc_.dt << "\n";
//         }

//         inline Path output_path() const
//         {
//             return desc_.output_path;
//         }

//     public:
//         std::shared_ptr<libMesh::DistributedMesh> mesh;
//         ProductFunctionSpace<LibMeshFunctionSpace> V;
//         std::shared_ptr<libMesh::EquationSystems> equation_systems;
//         std::shared_ptr<ElasticMaterial<USparseMatrix, UVector> > material;
//         std::shared_ptr<CompositeForcingFunction<UVector> > forcing_function;

//         int main_sys_num;
//         int aux_sys_num;

//         Desc desc_;
//     };

//     class ContactSimulation : public ElasticitySimulation {
//     public:
//         virtual ~ContactSimulation() {}

//         virtual bool init(libMesh::Parallel::Communicator &comm, Input &is) override
//         {
//             bool ok = true;
//             if(!ElasticitySimulation::init(comm, is)) {
//                 ok = false;
//             }

//             std::set<int> temp;
//             is.get("contact", [this,&temp](Input &is) {
//                 is.get("radius", contact_params.search_radius);

//                 std::string type;
//                 is.get("type", type);

//                 step_tol = 5e-6;
//                 is.get("step-tol", step_tol);

//                 max_nl_iter = 30;
//                 is.get("max-nl-iter", max_nl_iter);

//                 is_steady = false;
//                 n_transient_steps = 1;

//                 if(type == "steady") {
//                     is_steady = true;
//                 }

//                 use_pg = false;
//                 std::string solver;
//                 is.get("solver", solver);

//                 if(solver == "pg") {
//                     use_pg = true;
//                 }

//                 is.get("n-transient-steps", n_transient_steps);

//                 is.get("pairs", [this,&temp](Input &is) {
//                     is.get_all([this,&temp](Input &is) {
//                         int master = -1, slave = -1;
//                         is.get("master", master);
//                         is.get("slave", slave);

//                         // std::cout << master << " " << slave << std::endl;

//                         assert(master != -1);
//                         assert(slave  != -1);
//                         temp.insert(master);
//                         temp.insert(slave);

//                         contact_params.contact_pair_tags.push_back({ master, slave });
//                     });
//                 });
//             });

//             contact_surfaces.clear();
//             contact_surfaces.insert(contact_surfaces.end(), temp.begin(), temp.end());

//             return true;
//         }

//         virtual void describe(std::ostream &os) const override
//         {
//             ElasticitySimulation::describe(os);
//             contact_params.describe(os);
//         }

//         ContactParams contact_params;
//         std::vector<int> contact_surfaces;
//         bool is_steady;
//         int n_transient_steps;
//         double step_tol;
//         int max_nl_iter;
//         bool use_pg;

//     };

//     class WearSimulation::SimulationInput : public ContactSimulation {
//     public:
//         SimulationInput()
//         : wear_coefficient(7e-3), extrapolation_factor(10.)
//         {}

//         virtual bool init(libMesh::Parallel::Communicator &comm, Input &is) override
//         {
//             if(!ContactSimulation::init(comm, is)) return false;

//             is.get("wear", [this](Input &is) {
//                 is.get("n-cycles",   n_cycles);
//                 is.get("gait-cycle", gc);
//                 is.get("coeff", wear_coefficient);
//                 is.get("extrapolation", extrapolation_factor);
//             });

//             gc.init(V);
//             return true;
//         }

//         GaitCycle gc;
//         int n_cycles;
//         double wear_coefficient, extrapolation_factor;

//         virtual void describe(std::ostream &os) const override
//         {
//             ContactSimulation::describe(os);
//             contact_params.describe(os);

//             std::cout << "n_cycles: " << n_cycles << std::endl;
//             gc.describe(os);
//         }
//     };

//     void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &conf_file_path)
//     {
//         typedef utopia::ContactSolver<USparseMatrix, UVector> ContactSolverT;
//         typedef utopia::ContactStabilizedNewmark<USparseMatrix, UVector> TransientContactSolverT;

//         SimulationInput in;
//         auto is_ptr = open_istream(conf_file_path);
//         if(!is_ptr) {
//             std::cerr << "[Error] invalid path " << conf_file_path << std::endl;
//             assert(false);
//             return;
//         }

//         in.init_sim(init.comm(), *is_ptr);
//         in.describe(std::cout);

//         auto configuration_forces = std::make_shared<ConstantForcingFunction<UVector>>();
//         configuration_forces->value() = local_zeros(in.V.subspace(0).dof_map().n_local_dofs());
//         auto configuration_forced_material = std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(in.material,
//         configuration_forces);

//         std::shared_ptr<ContactSolverT> solver;
//         std::shared_ptr<TransientContactSolverT> transient_solver;
//         if(in.is_steady) {
//             solver = std::make_shared<ContactSolverT>(make_ref(in.V), configuration_forced_material,
//             in.contact_params);
//         } else {
//             double dt = in.gc.conf().dt() / in.n_transient_steps;
//             transient_solver = std::make_shared<TransientContactSolverT>(make_ref(in.V),
//             configuration_forced_material, dt, in.contact_params); solver = transient_solver;
//         }

//         std::cout << "n_dofs: " << in.V[0].dof_map().n_dofs() << std::endl;

//         // auto sor = std::make_shared<SOR<USparseMatrix, UVector>>();
//         // sor->rtol(1e-8);
//         // sor->atol(1e-16);
//         // sor->stol(1e-10);
//         // solver->set_linear_solver(sor);
//         solver->set_linear_solver(std::make_shared<Factorization<USparseMatrix, UVector>>());

//         solver->set_tol(in.step_tol);
//         solver->set_max_non_linear_iterations(in.max_nl_iter);
//         solver->set_max_outer_loops(40);

//         if(in.is_steady) {
//             if(in.use_pg) {
//                 solver->set_use_pg(true);
//             } else {
//                 solver->set_use_ssn(true);
//             }
//             // solver->tao().set_type("gpcg");
//             // solver->tao().verbose(true);
//             // solver->tao().atol(1e-8);
//             // solver->tao().rtol(1e-8);
//             // solver->tao().stol(1e-8);
//         } else {
//             // solver->tao().atol(1e-8);
//             // solver->tao().rtol(1e-8);
//             // solver->tao().stol(1e-8);
//             // solver->tao().verbose(true); //REMOVED_TRILINOS
//         }

//         // auto ls = std::make_shared<Factorization<USparseMatrix, UVector>>();
//         // auto ls = std::make_shared<GMRES<USparseMatrix, UVector>>();
//         // ls->atol(1e-15);
//         // ls->rtol(1e-15);
//         // ls->stol(1e-15);
//         // ls->max_it(1000);
//         // // ls->verbose(true);
//         // solver->set_linear_solver(ls);
//         // solver->set_bypass_contact(true);
//         // solver->set_max_outer_loops(10);

//         Wear wear;
//         wear.set_params(in.wear_coefficient, in.extrapolation_factor);
//         wear.init_aux_system(
//             *in.equation_systems
//         );

//         UVector overriden_displacement = local_zeros(in.V.subspace(0).dof_map().n_local_dofs());
//         UVector wear_displacement = overriden_displacement;

//         MechanicsState state;

//         libMesh::Nemesis_IO(*in.mesh).write_timestep(in.output_path() / "wear_in.e", *in.equation_systems, (1), (1));

//         libMesh::Nemesis_IO wear_io(*in.mesh);
//         libMesh::Nemesis_IO wear_ovv_io(*in.mesh);
//         libMesh::Nemesis_IO io(*in.mesh);

//         wear_io.write_timestep(in.output_path() / "wear.e", *in.equation_systems, 1, 0);

//         for(int i = 1; i <= in.n_cycles; ++i) {

//             std::cout << "cycle " << i << std::endl;

//             //gait-cycle
//             for(int t = 0; t < in.gc.conf().n_steps(); ++t) {
//                 std::cout << "\t step " << t << std::endl;
//                 configuration_forced_material->clear();
//                 in.gc.update(t);

//                 //set-up experiment
//                     //transform mesh

//                 overriden_displacement.set(0.);
//                 in.gc.displacement_and_forces(in.V, overriden_displacement, configuration_forces->value());

//                 {
//                     auto &sys = in.equation_systems->get_system<libMesh::LinearImplicitSystem>("wear");
//                     UVector temp = overriden_displacement + wear_displacement;
//                     convert(temp, *sys.solution);
//                     sys.solution->close();

//                     auto frame = ((i-1) * in.gc.conf().n_steps() + t + 1);
//                     wear_ovv_io.write_timestep(in.output_path() / "wear_overr.e", *in.equation_systems, frame,
//                     frame);
//                 }

//                 apply_displacement(overriden_displacement + wear_displacement, in.V.subspace(0).dof_map(), *in.mesh);

//                 //solve
//                 if(in.is_steady) {
//                     solver->solve_steady();
//                 } else {
//                     transient_solver->initial_condition(1.);
//                     transient_solver->solve_dynamic(in.n_transient_steps);
//                 }

//                 //update-wear
//                     //compute velocity
//                     //compute normal stress
//                     //compute sliding distance
//                     //compute wear

//                 state.displacement           = solver->displacement();
//                 state.displacement_increment = solver->displacement();

//                 state.velocity               = (1./in.gc.conf().dt()) * solver->displacement();

//                 if(in.forcing_function) {
//                     in.forcing_function->eval(state.displacement, state.external_force);
//                     state.external_force += configuration_forces->value();
//                 }

//                 solver->stress(state.displacement, state.stress);

//                 const double mag_stress   = norm2(state.stress);
//                 const double mag_velocity = norm2(state.velocity);

//                 // std::cout << "mag_stress:   " << mag_stress   << std::endl;
//                 // std::cout << "mag_velocity: " << mag_velocity << std::endl;

//                 //FIXME
//                 state.t = in.gc.conf().dt() * t;

//                 wear.update_aux_system(
//                     0,
//                     state,
//                     solver->contact(),
//                     in.gc.conf().dt(),
//                     *in.equation_systems
//                 );

//                 //clean-up experiment
//                     //transform-back mesh
//                 apply_displacement(-(overriden_displacement + wear_displacement), in.V.subspace(0).dof_map(),
//                 *in.mesh);

//                 auto &sys = in.equation_systems->get_system<libMesh::LinearImplicitSystem>("wear");
//                 UVector temp = state.displacement + overriden_displacement + wear_displacement;
//                 convert(temp, *sys.solution);
//                 sys.solution->close();

//                 auto frame = ((i-1) * in.gc.conf().n_steps() + t + 1);
//                 io.write_timestep(in.output_path() / "gait_cycles.e", *in.equation_systems, frame, frame);
//             }

//             //deform geometry
//             wear.mesh_displacement(in.V, in.contact_surfaces, wear_displacement);
//             wear_io.write_timestep(in.output_path() / "wear.e", *in.equation_systems, (i+1), i);
//         }

//         std::ofstream os(in.output_path() / "wear_profile.txt");

//         if(os.good()) {
//             wear.print(os);
//         }

//         os.close();

//         std::cout << "finished" << std::endl;
//     }

//     WearSimulation::WearSimulation()
//     {}

//     WearSimulation::~WearSimulation()
//     {}
// }

// #else

// namespace utopia {
//     void WearSimulation::run(libMesh::LibMeshInit &init, const std::string &conf_file_path)
//     {
//       std::cerr << "[Error] not run (does not work with trilinos backend)" << std::endl;
//     }

//     WearSimulation::WearSimulation()
//     {}

//     WearSimulation::~WearSimulation()
//     {}
// }

// #endif //UTOPIA_WITH_TRILINOS_ALGEBRA
