#include "utopia_FractureFlowApp.hpp"

#include "utopia_TransferApp.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_MeshTransferOperator.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_Blocks.hpp"
#include "utopia_Eval_Blocks.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"

namespace utopia {
    
    
    typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
    
    void FractureFlowApp::init(libMesh::LibMeshInit &init)
    {
        comm_ = make_ref(init.comm());
    }


    static void refine_around_fractures(
        const std::shared_ptr<libMesh::UnstructuredMesh> &fracture_network,
        const libMesh::Order &elem_order,
        const std::shared_ptr<libMesh::UnstructuredMesh> &mesh,
        const int refinement_loops = 1,
        const bool use_interpolation = false
        )
    {

        libMesh::MeshRefinement mesh_refinement(*mesh);

        for(int i = 0; i < refinement_loops; ++i) {
        //equations system
            auto vol_equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
            auto &vol_sys = vol_equation_systems->add_system<libMesh::LinearImplicitSystem>("vol_sys");

            auto surf_equation_systems = std::make_shared<libMesh::EquationSystems>(*fracture_network);
            auto &surf_sys = surf_equation_systems->add_system<libMesh::LinearImplicitSystem>("surf_sys");

        //scalar function space
            auto V_vol  = FunctionSpaceT(vol_equation_systems, libMesh::LAGRANGE, elem_order,      "u_vol");
            auto V_surf = FunctionSpaceT(surf_equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_surf");

            V_vol.initialize();
            V_surf.initialize();

            Chrono c;
            c.start();
            USparseMatrix B;
            moonolith::Communicator comm(mesh->comm().get());
            if(assemble_volume_transfer(
                comm,
                mesh,
                fracture_network,
                make_ref(V_vol.dof_map()),
                make_ref(V_surf.dof_map()),
                0,
                0,
                true,
                1,
                B,
                {},
                use_interpolation))
            {
                c.stop();
                std::cout << c << std::endl;

                Interpolator interp(make_ref(B));
                interp.normalize_rows();
                interp.describe(std::cout);

                USparseMatrix D_inv = diag(1./sum(B, 1));
                USparseMatrix T = D_inv * B;

                USparseMatrix T_t = transpose(T);
                UVector t_temp = sum(T_t, 1);
                UVector t = ghosted(local_size(t_temp).get(0), size(t_temp).get(0), V_vol.dof_map().get_send_list());
                t = t_temp;


                std::vector<libMesh::dof_id_type> indices;
                std::vector<double> values;

                mesh_refinement.clean_refinement_flags();

                Read<UVector> r_(t);
                for(auto e_it = elements_begin(*mesh); e_it != elements_end(*mesh); ++e_it) {
                    V_vol.dof_map().dof_indices(*e_it, indices);
                    t.get(indices, values);

                    double val = std::accumulate(values.begin(), values.end(), 0.,  std::plus<double>());
                    if(val > 0) {
                        (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
                    }

                }

                mesh_refinement.make_flags_parallel_consistent();
                mesh_refinement.refine_elements();
                mesh_refinement.test_level_one(true);

            } else {
                assert(false);
            }
        }

        // mesh_refinement.clean_refinement_flags();
        mesh->prepare_for_use();
    }

    static bool assemble_interpolation(FunctionSpaceT &from, FunctionSpaceT &to, USparseMatrix &B, USparseMatrix &D)
    {
        auto assembler = std::make_shared<InterpolationLocalAssembler>(from.mesh().mesh_dimension());
        auto local2global = std::make_shared<Local2Global>(true);
        
        TransferAssembler transfer_assembler(assembler, local2global);
        
        std::vector< std::shared_ptr<USparseMatrix> > mats;
        if(!transfer_assembler.assemble(
                                        make_ref(from.mesh()),
                                        make_ref(from.dof_map()),
                                        make_ref(to.mesh()),
                                        make_ref(to.dof_map()),
                                        mats)) {
            return false;
        }
        
        B = std::move(*mats[0]);
        D = diag(sum(B, 1));

        double sum_B = sum(B);
        double sum_D = sum(D);

        std::cout << sum_B << " == " << sum_D << std::endl;
        return true;
    }
    
    
    static bool assemble_projection(FunctionSpaceT &from, FunctionSpaceT &to, USparseMatrix &B, USparseMatrix &D)
    {
        auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), false, true);
        auto local2global = std::make_shared<Local2Global>(false);
        
        TransferAssembler transfer_assembler(assembler, local2global);
        
        std::vector< std::shared_ptr<USparseMatrix> > mats;
        if(!transfer_assembler.assemble(
                                        make_ref(from.mesh()),
                                        make_ref(from.dof_map()),
                                        make_ref(to.mesh()),
                                        make_ref(to.dof_map()),
                                        mats)) {
            return false;
        }
        
        B = std::move(*mats[0]);
        D = std::move(*mats[1]);

        double sum_B = sum(B);
        double sum_D = sum(D);

        std::cout << sum_B << " == " << sum_D << std::endl;
        return true;
    }


    static bool assemble_coupling(FunctionSpaceT &from, FunctionSpaceT &to, USparseMatrix &B)
    {
        auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), false, false);
        auto local2global = std::make_shared<Local2Global>(false);
        
        TransferAssembler transfer_assembler(assembler, local2global);
        
        std::vector< std::shared_ptr<USparseMatrix> > mats;
        if(!transfer_assembler.assemble(
                                        make_ref(from.mesh()),
                                        make_ref(from.dof_map()),
                                        make_ref(to.mesh()),
                                        make_ref(to.dof_map()),
                                        mats)) {
            return false;
        }
        
        B = std::move(*mats[0]);

        double sum_B = sum(B);

        std::cout << sum_B << std::endl;
        return true;
    }

    //use different lagr mult space
    static bool assemble_projection(
        FunctionSpaceT &from,
        FunctionSpaceT &to, 
        FunctionSpaceT &lagr, 
        USparseMatrix &B, USparseMatrix &D)
    {
        if(assemble_coupling(from, lagr, B)) {
            return assemble_coupling(to, lagr, D);
        } else {
            return false;
        }
    }
    
    static void solve_monolithic(FunctionSpaceT &V_m,
                                 FunctionSpaceT &V_s,
                                 USparseMatrix &A_m,
                                 UVector &rhs_m,
                                 USparseMatrix &A_s,
                                 UVector &rhs_s,
                                 UVector &sol_m,
                                 UVector &sol_s,
                                 UVector &lagr)
    {
        
        USparseMatrix B, D;
        assemble_projection(V_m, V_s, B, D);
        // assemble_interpolation(V_m, V_s, B, D);

        D *= -1.;

        USparseMatrix D_t = transpose(D);
        USparseMatrix B_t = transpose(B);

        set_zero_at_constraint_rows(V_m.dof_map(), B_t);
        set_zero_at_constraint_rows(V_s.dof_map(), D_t);

        USparseMatrix A = Blocks<USparseMatrix>(3, 3, 
        {
            make_ref(A_m), nullptr, make_ref(B_t),
            nullptr, make_ref(A_s), make_ref(D_t),
            make_ref(B), make_ref(D), nullptr
        });

        UVector z = local_zeros(local_size(rhs_s));
        UVector rhs = blocks(rhs_m, rhs_s, z);

        sol_m = local_zeros(local_size(rhs_m));
        sol_s = local_zeros(local_size(rhs_s));
        lagr  = local_zeros(local_size(rhs_s));

        UVector sol = blocks(sol_m, sol_s, lagr);

        Factorization<USparseMatrix, UVector> op;
        op.update(make_ref(A));
        op.apply(rhs, sol);

        //write("A.m", A);

        undo_blocks(sol, sol_m, sol_s, lagr);
    }


    static void solve_monolithic(FunctionSpaceT &V_m,
                                 FunctionSpaceT &V_s,
                                 FunctionSpaceT &L,
                                 USparseMatrix &A_m,
                                 UVector &rhs_m,
                                 USparseMatrix &A_s,
                                 UVector &rhs_s,
                                 UVector &sol_m,
                                 UVector &sol_s,
                                 UVector &lagr)
    {
        
        USparseMatrix B, D;
        assemble_projection(V_m, V_s, L, B, D);
        // assemble_interpolation(V_m, V_s, B, D);

        D *= -1.;

        USparseMatrix D_t = transpose(D);
        USparseMatrix B_t = transpose(B);

        set_zero_at_constraint_rows(V_m.dof_map(), B_t);
        set_zero_at_constraint_rows(V_s.dof_map(), D_t);

        USparseMatrix A = Blocks<USparseMatrix>(3, 3, 
        {
            make_ref(A_m), nullptr, make_ref(B_t),
            nullptr, make_ref(A_s), make_ref(D_t),
            make_ref(B), make_ref(D), nullptr
        });

        UVector z = local_zeros(L.dof_map().n_local_dofs());
        UVector rhs = blocks(rhs_m, rhs_s, z);

        sol_m = local_zeros(local_size(rhs_m));
        sol_s = local_zeros(local_size(rhs_s));
        lagr  = local_zeros(L.dof_map().n_local_dofs());

        UVector sol = blocks(sol_m, sol_s, lagr);

        Factorization<USparseMatrix, UVector> op;
        op.solve(A, rhs, sol);

        undo_blocks(sol, sol_m, sol_s, lagr);
    }
    
    static void write_solution(const std::string &name,
                               UVector &sol,
                               const int time_step,
                               const double t,
                               FunctionSpaceT &space,
                               libMesh::Nemesis_IO &io)
    {
        utopia::convert(sol, *space.equation_system().solution);
        space.equation_system().solution->close();
        io.write_timestep(name, space.equation_systems(), time_step, t);
    }
    
    static void solve_staggered(const std::string &operator_type,
                                FunctionSpaceT &V_m,
                                FunctionSpaceT &V_s,
                                USparseMatrix &A_m,
                                UVector &rhs_m,
                                USparseMatrix &A_s,
                                UVector &rhs_s,
                                UVector &sol_m,
                                UVector &sol_s,
                                UVector &lagr)
    {
        
        libMesh::Nemesis_IO io_m(V_m.mesh());
        libMesh::Nemesis_IO io_s(V_s.mesh());
        libMesh::Nemesis_IO io_m_l(V_m.mesh());
        libMesh::Nemesis_IO io_s_l(V_s.mesh());
        
        Factorization<USparseMatrix, UVector> op_m;
        op_m.update(make_ref(A_m));
        
        if(empty(sol_s)) {
            sol_s = local_zeros(local_size(rhs_s));
        }
        
        UVector lagr_m = local_zeros(local_size(rhs_m));
        UVector lagr_s = local_zeros(local_size(rhs_s));
        
        UVector rhs_lagr_m, rhs_lagr_s;
        UVector delta_lagr = local_zeros(local_size(rhs_s));
        
        double dumping = 1.;
        
        TransferOptions opts;
        opts.from_var_num = 0;
        opts.to_var_num   = 0;
        
        MeshTransferOperator t(make_ref(V_m.mesh()),
                               make_ref(V_m.dof_map()),
                               make_ref(V_s.mesh()),
                               make_ref(V_s.dof_map()),
                               opts
                               );
        
        t.initialize(operator_type);
        t.write("./");

        lagr_m.set(0.);
        
        for(int i = 0; i < 20; ++i) {
            apply_zero_boundary_conditions(V_m.dof_map(), lagr_m);
            rhs_lagr_m = rhs_m + lagr_m;
            
            write_solution("lagr_m.e",
                           lagr_m,
                           (i + 1),
                           (i + 1),
                           V_m,
                           io_m_l
                           );
            
            op_m.apply(rhs_lagr_m, sol_m);
            
            write_solution("sol_m.e",
                           sol_m,
                           (i + 1),
                           (i + 1),
                           V_m,
                           io_m
                           );
            
            sol_s.set(0.);
            t.apply(sol_m, sol_s);
            
            apply_boundary_conditions(V_s, A_s, sol_s);
            lagr_s = rhs_s - A_s * sol_s;
            apply_zero_boundary_conditions(V_s.dof_map(), lagr_s);

            
            write_solution("sol_s.e",
                           sol_s,
                           (i + 1),
                           (i + 1),
                           V_s,
                           io_s
                           );
            
            
            write_solution("lagr_s.e",
                           lagr_s,
                           (i + 1),
                           (i + 1),
                           V_s,
                           io_s_l
                           );
            
            double n_lagr_s = norm2(lagr_s);
            
            disp(n_lagr_s);
            
            if(n_lagr_s < 1e-14) {
                break;
            }
            
            delta_lagr.set(0);
            t.apply_transpose(lagr_s, delta_lagr);
            lagr_m += dumping * delta_lagr;
        }
        
        lagr = lagr_m;
    }
    
    static void solve_separate(const std::string &operator_type,
                               FunctionSpaceT &V_m,
                               FunctionSpaceT &V_s,
                               USparseMatrix &A_m,
                               UVector &rhs_m,
                               USparseMatrix &A_s,
                               UVector &rhs_s,
                               UVector &sol_m,
                               UVector &sol_s,
                               UVector &lagr)
    {
        Factorization<USparseMatrix, UVector> op_m;
        op_m.update(make_ref(A_m));
        
        Factorization<USparseMatrix, UVector> op_s;
        op_s.update(make_ref(A_s));
        
        
        op_m.apply(rhs_m, sol_m);
        op_s.apply(rhs_s, sol_s);
        
        
        MeshTransferOperator t(make_ref(V_m.mesh()),
                               make_ref(V_m.dof_map()),
                               make_ref(V_s.mesh()),
                               make_ref(V_s.dof_map())
                               );
        
        t.initialize(operator_type);
        
        UVector sol_transfered;
        t.apply(sol_m, sol_transfered);
        lagr = sol_transfered - sol_s;
    }
    
    
    static void refine(const int n_refs, libMesh::MeshBase &mesh)
    {
        if(n_refs <= 0) return;
        
        libMesh::MeshRefinement mesh_refinement(mesh);
        mesh_refinement.make_flags_parallel_consistent();
        mesh_refinement.uniformly_refine(n_refs);
    }
    
    class FractureFlowApp::Input : public Serializable {
    public:
        Input(libMesh::Parallel::Communicator &comm)
        : mesh(comm), space(mesh)
        {}

        void read(InputStream &is) override
        {
            try {
                is.read("mesh", mesh);
                is.read("space", space);

                forcing_function = std::make_shared< UIForcingFunction<FunctionSpaceT, UVector> >(space.subspace(0));
                is.read("forcing-function", *forcing_function);
    
                //material parameters              
                is.read("diffusivity", diffusivity);

            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
            }
        }

        inline bool empty() const
        {
            return mesh.empty();
        }

        void describe(std::ostream &os = std::cout) const
        {
            // os << "-----------------------------------\n";
            // mesh.describe(os);
            // space.describe(os);
            // forcing_function.describe(os);
            // os << "-----------------------------------\n";
        }
        
        UIMesh<libMesh::DistributedMesh> mesh;    
        UIFunctionSpace<FunctionSpaceT>  space;
        std::shared_ptr< UIForcingFunction<FunctionSpaceT, UVector> > forcing_function;

        double diffusivity;
    };
    
    void FractureFlowApp::run(const std::string &conf_file_path)
    {
        Chrono c;
        
        c.start();


        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////// SET-UP ////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        auto is_ptr = open_istream(conf_file_path);
        
        Input master_in(*comm_), slave_in(*comm_);
        // Input multiplier_in(*comm_); //LAMBDA
        
        is_ptr->read("master", master_in);
        is_ptr->read("slave",  slave_in);

        // is_ptr->read("multiplier", multiplier_in); //LAMBDA

        // if(multiplier_in.empty()) { //LAMBDA
            // multiplier_in = slave_in;  //LAMBDA
        // }                           //LAMBDA
        
        std::string solve_strategy = "monolithic";
        is_ptr->read("solve-strategy", solve_strategy);


        std::string operator_type = "L2_PROJECTION";
        is_ptr->read("operator-type", operator_type);
        
        master_in.describe();
        slave_in.describe();
        // multiplier_in.describe(); //LAMBDA
        
        std::cout << "solve_strategy: "  << solve_strategy << std::endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
            
        auto &V_m = master_in.space.subspace(0);
        auto &V_s = slave_in.space.subspace(0);

        //////////////////////////// Variational formulation ////////////////////////////
        
        auto u_m = trial(V_m);
        auto v_m = test(V_m);
        
        auto u_s = trial(V_s);
        auto v_s = test(V_s);
        
        // L.initialize(); //LAMBDA

        std::cout << "n_dofs: " << V_m.dof_map().n_dofs() << " x " <<  V_s.dof_map().n_dofs();
        // std::cout << " x " <<  L.dof_map().n_dofs() << ""; //LAMBDA
        std::cout << std::endl;
        
        auto eq_m = master_in.diffusivity * inner(grad(u_m), grad(v_m)) * dX;
        auto eq_s = slave_in.diffusivity  * inner(grad(u_s), grad(v_s)) * dX;
        
        //////////////////////////// Generation of the algebraic system ////////////////////////////

        USparseMatrix A_m, A_s;
        UVector rhs_m, rhs_s;
        utopia::assemble(eq_m, A_m);
        utopia::assemble(eq_s, A_s);

        UVector x_m = local_zeros(V_m.dof_map().n_local_dofs());
        UVector x_s = local_zeros(V_s.dof_map().n_local_dofs());

        master_in.forcing_function->eval(x_m, rhs_m);
        slave_in.forcing_function->eval(x_s, rhs_s);
        
        apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
        apply_boundary_conditions(V_s.dof_map(), A_s, rhs_s);
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
       
        UVector sol_m, sol_s, lagr;

        if(solve_strategy == "staggered") {
            
            solve_staggered(operator_type,
                            V_m,
                            V_s,
                            A_m,
                            rhs_m,
                            A_s,
                            rhs_s,
                            sol_m,
                            sol_s,
                            lagr
                            );
            
        } else {
            solve_monolithic(V_m,
                             V_s,
                             // L, //LAMBDA
                             A_m,
                             rhs_m,
                             A_s,
                             rhs_s,
                             sol_m,
                             sol_s,
                             lagr
                             );
        }
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        

        // UGLY code for writing stuff to disk
        libMesh::Nemesis_IO io_m(master_in.mesh.mesh());
        libMesh::Nemesis_IO io_s(slave_in.mesh.mesh());
        // libMesh::Nemesis_IO io_multiplier(*mesh_multiplier);
        
        utopia::convert(sol_m, *V_m.equation_system().solution);
        V_m.equation_system().solution->close();
        io_m.write_timestep(V_m.equation_system().name() + ".e", V_m.equation_systems(), 1, 0);
        
        
        utopia::convert(sol_s, *V_s.equation_system().solution);
        V_s.equation_system().solution->close();
        io_s.write_timestep(V_s.equation_system().name() + ".e", V_s.equation_systems(), 1, 0);

        // utopia::convert(lagr, *L.equation_system().solution);                                         //LAMBDA
        // L.equation_system().solution->close();                                                        //LAMBDA
        // io_multiplier.write_timestep(L.equation_system().name() + ".e", L.equation_systems(), 1, 0);  //LAMBDA
        
        
        c.stop();
        std::cout << c << std::endl;
    }
}

