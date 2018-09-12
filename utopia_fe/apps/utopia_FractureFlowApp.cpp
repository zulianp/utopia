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

#include "libmesh/mesh_refinement.h"

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
            DSMatrixd B;
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

                DSMatrixd D_inv = diag(1./sum(B, 1));
                DSMatrixd T = D_inv * B;

                DSMatrixd T_t = transpose(T);
                DVectord t_temp = sum(T_t, 1);
                DVectord t = ghosted(local_size(t_temp).get(0), size(t_temp).get(0), V_vol.dof_map().get_send_list());
                t = t_temp;


                std::vector<libMesh::dof_id_type> indices;
                std::vector<double> values;

                mesh_refinement.clean_refinement_flags();

                Read<DVectord> r_(t);
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

    static bool assemble_interpolation(FunctionSpaceT &from, FunctionSpaceT &to, DSMatrixd &B, DSMatrixd &D)
    {
        auto assembler = std::make_shared<InterpolationLocalAssembler>(from.mesh().mesh_dimension());
        auto local2global = std::make_shared<Local2Global>(true);
        
        TransferAssembler transfer_assembler(assembler, local2global);
        
        std::vector< std::shared_ptr<DSMatrixd> > mats;
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
    
    
    static bool assemble_projection(FunctionSpaceT &from, FunctionSpaceT &to, DSMatrixd &B, DSMatrixd &D)
    {
        auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), false, true);
        auto local2global = std::make_shared<Local2Global>(false);
        
        TransferAssembler transfer_assembler(assembler, local2global);
        
        std::vector< std::shared_ptr<DSMatrixd> > mats;
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

    //use different lagr mult space
    static bool assemble_projection(
        FunctionSpaceT &from,
        FunctionSpaceT &to, 
        FunctionSpaceT &lagr, 
        DSMatrixd &B, DSMatrixd &D)
    {

        {
            DSMatrixd dump;
            assemble_projection(from, lagr, B, dump);
        }

        {
            DSMatrixd dump;
            assemble_projection(to, lagr, D, dump);
        }

        return true;
    }
    
    
    // static void solve_monolithic_old(FunctionSpaceT &V_m,
    //                              FunctionSpaceT &V_s,
    //                              DSMatrixd &A_m,
    //                              DVectord &rhs_m,
    //                              DSMatrixd &A_s,
    //                              DVectord &rhs_s,
    //                              DVectord &sol_m,
    //                              DVectord &sol_s,
    //                              DVectord &lagr)
    // {
        
    //     DSMatrixd B, D;
    //     assemble_projection(V_m, V_s, B, D);
    //     // assemble_interpolation(V_m, V_s, B, D);

    //     D *= -1.;
    //     // B *= -1.;
        
    //     auto s_m = local_size(A_m);
    //     auto s_s = local_size(A_s);
        
        
    //     auto nnz_x_row_m =
    //     std::max(*std::max_element(V_m.dof_map().get_n_nz().begin(), V_m.dof_map().get_n_nz().end()),
    //              *std::max_element(V_m.dof_map().get_n_oz().begin(), V_m.dof_map().get_n_oz().end()));
        
    //     auto nnz_x_row_s =
    //     std::max(*std::max_element(V_s.dof_map().get_n_nz().begin(), V_s.dof_map().get_n_nz().end()),
    //              *std::max_element(V_s.dof_map().get_n_oz().begin(), V_s.dof_map().get_n_oz().end()));
        
        
    //     DSMatrixd A = local_sparse(
    //                                s_m.get(0) + 2 * s_s.get(0),
    //                                s_m.get(1) + 2 * s_s.get(1),
    //                                2 * nnz_x_row_m + 2 * nnz_x_row_s //FIXME
    //                                );
        
        
    //     auto rr_m = row_range(A_m);
    //     auto rr_s = row_range(A_s);
        
    //     auto off_r   = size(A_m).get(0);
    //     auto off_c   = size(A_m).get(1);
        
    //     auto off_r_l = off_r + size(A_s).get(0);
    //     auto off_c_l = off_c + size(A_s).get(1);
        
    //     {
    //         Write<DSMatrixd> w_(A);
            
    //         for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
    //             {
    //                 RowView<DSMatrixd> row(A_m, i);
                    
    //                 for(auto k = 0; k < row.n_values(); ++k) {
    //                     //(1, 1)
    //                     A.set(i, row.col(k), row.get(k));
    //                 }
    //             }
    //         }
            
    //         for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
                
    //             {
    //                 RowView<DSMatrixd> row(A_s, i);
                    
    //                 for(auto k = 0; k < row.n_values(); ++k) {
    //                     //(2, 2)
    //                     A.set(off_r + i, off_c + row.col(k), row.get(k));
    //                 }
    //             }
                
    //             {
    //                 RowView<DSMatrixd> row(B, i);
                    
    //                 for(auto k = 0; k < row.n_values(); ++k) {
    //                     //(1, 3)  //B^T
    //                     if(!V_m.dof_map().is_constrained_dof(row.col(k))) {
    //                         A.set(row.col(k), off_c_l + i, row.get(k));
    //                     }
                        
    //                     //(3, 1) //B
    //                     A.set(off_r_l + i, row.col(k), row.get(k));
    //                 }
    //             }
                
    //             {
    //                 RowView<DSMatrixd> row(D, i);
                    
    //                 for(auto k = 0; k < row.n_values(); ++k) {
                        
    //                     //(2, 3) //D^T
    //                     if(!V_s.dof_map().is_constrained_dof(row.col(k))) {
    //                         A.set(off_r + row.col(k), off_c_l + i, row.get(k));
    //                     }
                        
    //                     //(3, 2)   //D
    //                     A.set(off_r_l + i, off_c + row.col(k), row.get(k));
    //                 }
    //             }
                
    //         }
    //     }
        
    //     DVectord rhs = local_zeros(local_size(rhs_m).get(0) + 2 * local_size(rhs_s).get(0));
        
    //     {
    //         Write<DVectord> w_(rhs);
    //         Write<DVectord> r_m(rhs_m), r_s(rhs_s);
            
    //         for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
    //             rhs.set(i, rhs_m.get(i));
    //         }
            
    //         for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
    //             rhs.set(off_r + i, rhs_s.get(i));
    //         }
    //     }
        
    //     // A.implementation().set_name("a");
    //     // write("A.m", A);
        
    //     Factorization<DSMatrixd, DVectord> op;
    //     // GMRES<DSMatrixd, DVectord> op;
    //     // LUDecomposition<DSMatrixd, DVectord> op;
    //     op.update(make_ref(A));
        
        
    //     DVectord sol = local_zeros(local_size(rhs));
    //     op.apply(rhs, sol);
        
        
    //     sol_m  = local_zeros(local_size(rhs_m));
    //     sol_s  = local_zeros(local_size(rhs_s));
    //     lagr   = local_zeros(local_size(rhs_s));
        
    //     {
    //         Write<DVectord> w_(sol);
    //         Write<DVectord> r_m(sol_m), r_s(sol_s);
            
    //         for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
    //             sol_m.set(i, sol.get(i));
    //         }
            
    //         for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
    //             sol_s.set(i, sol.get(off_r + i));
    //         }
            
    //         for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
    //             lagr.set(i, sol.get(off_r_l + i));
    //         }
    //     }
        
    // }


    static void solve_monolithic(FunctionSpaceT &V_m,
                                 FunctionSpaceT &V_s,
                                 DSMatrixd &A_m,
                                 DVectord &rhs_m,
                                 DSMatrixd &A_s,
                                 DVectord &rhs_s,
                                 DVectord &sol_m,
                                 DVectord &sol_s,
                                 DVectord &lagr)
    {
        
        DSMatrixd B, D;
        assemble_projection(V_m, V_s, B, D);
        // assemble_interpolation(V_m, V_s, B, D);

        D *= -1.;

        DSMatrixd D_t = transpose(D);
        DSMatrixd B_t = transpose(B);

        set_zero_at_constraint_rows(V_m.dof_map(), B_t);
        set_zero_at_constraint_rows(V_s.dof_map(), D_t);

        DSMatrixd A = Blocks<DSMatrixd>(3, 3, 
        {
            make_ref(A_m), nullptr, make_ref(B_t),
            nullptr, make_ref(A_s), make_ref(D_t),
            make_ref(B), make_ref(D), nullptr
        });

        DVectord z = local_zeros(local_size(rhs_s));

        DVectord rhs = Blocks<DVectord>(
        {
            make_ref(rhs_m),
            make_ref(rhs_s),
            make_ref(z)
        });

        sol_m = local_zeros(local_size(rhs_m));
        sol_s = local_zeros(local_size(rhs_s));
        lagr  = local_zeros(local_size(rhs_s));

        DVectord sol = Blocks<DVectord>(
        {
            make_ref(sol_m),
            make_ref(sol_s),
            make_ref(lagr)
        });

        Factorization<DSMatrixd, DVectord> op;
        op.update(make_ref(A));
        op.apply(rhs, sol);

        auto r_m = range(sol_m);
        auto r_s = range(sol_s);

        auto r = range(sol);

        auto n_m = local_size(sol_m).get(0);
        auto n_s = n_m + local_size(sol_s).get(0);

        {
            Read<DVectord> r_(sol);
            Write<DVectord> w_m(sol_m), w_s(sol_s), w_l(lagr);

            SizeType index = r.begin();
            for(auto i = r_m.begin(); i < r_m.end(); ++i) {
                sol_m.set(i, sol.get(index++));
            }

            for(auto i = r_s.begin(); i < r_s.end(); ++i) {
                sol_s.set(i, sol.get(index++));
            }

            for(auto i = r_s.begin(); i < r_s.end(); ++i) {
                lagr.set(i, sol.get(index++));
            }
        }
    }
    
    static void write_solution(const std::string &name,
                               DVectord &sol,
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
                                DSMatrixd &A_m,
                                DVectord &rhs_m,
                                DSMatrixd &A_s,
                                DVectord &rhs_s,
                                DVectord &sol_m,
                                DVectord &sol_s,
                                DVectord &lagr)
    {
        
        libMesh::Nemesis_IO io_m(V_m.mesh());
        libMesh::Nemesis_IO io_s(V_s.mesh());
        libMesh::Nemesis_IO io_m_l(V_m.mesh());
        libMesh::Nemesis_IO io_s_l(V_s.mesh());
        
        Factorization<DSMatrixd, DVectord> op_m;
        op_m.update(make_ref(A_m));
        
        if(empty(sol_s)) {
            sol_s = local_zeros(local_size(rhs_s));
        }
        
        DVectord lagr_m = local_zeros(local_size(rhs_m));
        DVectord lagr_s = local_zeros(local_size(rhs_s));
        
        DVectord rhs_lagr_m, rhs_lagr_s;
        DVectord delta_lagr = local_zeros(local_size(rhs_s));
        
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
                               DSMatrixd &A_m,
                               DVectord &rhs_m,
                               DSMatrixd &A_s,
                               DVectord &rhs_s,
                               DVectord &sol_m,
                               DVectord &sol_s,
                               DVectord &lagr)
    {
        Factorization<DSMatrixd, DVectord> op_m;
        op_m.update(make_ref(A_m));
        
        Factorization<DSMatrixd, DVectord> op_s;
        op_s.update(make_ref(A_s));
        
        
        op_m.apply(rhs_m, sol_m);
        op_s.apply(rhs_s, sol_s);
        
        
        MeshTransferOperator t(make_ref(V_m.mesh()),
                               make_ref(V_m.dof_map()),
                               make_ref(V_s.mesh()),
                               make_ref(V_s.dof_map())
                               );
        
        t.initialize(operator_type);
        
        DVectord sol_transfered;
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
        
        void read(InputStream &is) override
        {
            try {
                is.read("mesh", mesh_type);
                if(mesh_type == "file") {
                    is.read("path", path);
                }
                
                is.read("boundary-conditions", [this](InputStream &is) {
                    is.read_all([this](InputStream &is) {
                        int side_set = 0;
                        
                        is.read("side", side_set);
                        
                        
                        double value = 0;
                        is.read("value", value);
                        
                        sides.push_back(side_set);
                        values.push_back(value);
                    });
                });
                
                diffusivity = 1.;
                is.read("diffusivity", diffusivity);
                forcing_function = 0.;
                
                is.read("forcing-function", forcing_function);

                refinements = 0;
                is.read("refinements", refinements);


                adaptive_refinements = 0;
                is.read("adaptive-refinements", adaptive_refinements);

            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
            }
        }
        
        void make_mesh(libMesh::DistributedMesh &mesh) const
        {
            if(this->mesh_type == "file") {
                mesh.read(path);
            } else if(this->mesh_type == "unit-square") {
                const auto elem_type = libMesh::QUAD8;
                libMesh::MeshTools::Generation::build_square(mesh,
                                                             5, 5,
                                                             -0., 1.,
                                                             -0., 1.,
                                                             elem_type
                                                             );
            }
            
            refine(refinements, mesh);

        }
        
        void set_up_bc(const FunctionSpaceT &V)
        {
            auto u = trial(V);
            std::size_t n = sides.size();
            
            for(std::size_t i = 0; i < n; ++i) {
                init_constraints(
                                 constraints(
                                             boundary_conditions(u == coeff(values[i]), {sides[i]})
                                             )
                                 );
            }
        }
        
        void apply_adaptive_refinement(
            const std::shared_ptr<libMesh::UnstructuredMesh> &fracture_network,
            const libMesh::Order &elem_order,
            const std::shared_ptr<libMesh::UnstructuredMesh> &mesh)
        {
            if(adaptive_refinements) {
                refine_around_fractures(fracture_network, elem_order, mesh, adaptive_refinements, false);
            }
        }

        void desribe(std::ostream &os = std::cout) const
        {
            os << "-----------------------------------\n";
            os << "mesh_type:       " << mesh_type << "\n";
            os << "path:            " << path << "\n";
            os << "diffusivity:     " << diffusivity << "\n";
            os << "forcing_function: " << forcing_function << "\n";
            
            os << "side, value\n";
            
            std::size_t n = sides.size();
            for(std::size_t i = 0; i < n; ++i) {
                os << sides[i]<< ", " << values[i] << "\n";
            }
            
            os << "-----------------------------------\n";
        }
        
        
        std::string mesh_type, path;
        std::vector<int> sides;
        std::vector<double> values;
        double diffusivity;
        double forcing_function;
        int refinements;
        int adaptive_refinements;
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
        
        Input master_in, slave_in;
        
        is_ptr->read("master", master_in);
        is_ptr->read("slave", slave_in);
        
        std::string solve_strategy = "staggered";
        is_ptr->read("solve-strategy", solve_strategy);


        std::string operator_type = "INTERPOLATION";
        is_ptr->read("operator-type", operator_type);
        
        master_in.desribe();
        slave_in.desribe();
        
        std::cout << "solve_strategy: "  << solve_strategy << std::endl;
       
        auto mesh_master = std::make_shared<libMesh::DistributedMesh>(*comm_);
        auto mesh_slave = std::make_shared<libMesh::DistributedMesh>(*comm_);
        
        master_in.make_mesh(*mesh_master);
        slave_in.make_mesh(*mesh_slave);

        master_in.apply_adaptive_refinement(mesh_slave, libMesh::FIRST, mesh_master);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        //////////////////////////// Variational formulation ////////////////////////////
        //equations system
        auto equation_systems_master = std::make_shared<libMesh::EquationSystems>(*mesh_master);
        auto &sys_master = equation_systems_master->add_system<libMesh::LinearImplicitSystem>("master");
        
        auto equation_systems_slave = std::make_shared<libMesh::EquationSystems>(*mesh_slave);
        auto &sys_slave = equation_systems_slave->add_system<libMesh::LinearImplicitSystem>("slave");
        
        //scalar function space
        const auto elem_order = libMesh::FIRST;
        auto V_m = FunctionSpaceT(equation_systems_master, libMesh::LAGRANGE, elem_order, "u_m");
        auto V_s = FunctionSpaceT(equation_systems_slave, libMesh::LAGRANGE, elem_order,  "u_s");
        
        auto u_m = trial(V_m);
        auto v_m = test(V_m);
        
        auto u_s = trial(V_s);
        auto v_s = test(V_s);
        
        master_in.set_up_bc(V_m);
        slave_in.set_up_bc(V_s);
        
        V_m.initialize();
        V_s.initialize();
        
        auto eq_m = master_in.diffusivity * inner(grad(u_m), grad(v_m)) * dX == inner(coeff(master_in.forcing_function), v_m) * dX;
        auto eq_s = slave_in.diffusivity  * inner(grad(u_s), grad(v_s)) * dX == inner(coeff(slave_in.forcing_function),  v_s) * dX;
        
        //////////////////////////// Generation of the algebraic system ////////////////////////////

        DSMatrixd A_m, A_s;
        DVectord rhs_m, rhs_s;
        utopia::assemble(eq_m, A_m, rhs_m);
        utopia::assemble(eq_s, A_s, rhs_s);
        
        apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
        apply_boundary_conditions(V_s.dof_map(), A_s, rhs_s);
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
       
        DVectord sol_m, sol_s, lagr;

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
        
        libMesh::Nemesis_IO io_m(*mesh_master);
        libMesh::Nemesis_IO io_s(*mesh_slave);
        
        
        utopia::convert(sol_m, *V_m.equation_system().solution);
        V_m.equation_system().solution->close();
        io_m.write_timestep(V_m.equation_system().name() + ".e", V_m.equation_systems(), 1, 0);
        
        
        utopia::convert(sol_s, *V_s.equation_system().solution);
        V_s.equation_system().solution->close();
        io_s.write_timestep(V_s.equation_system().name() + ".e", V_s.equation_systems(), 1, 0);
        
        if(size(lagr) == size(sol_m)) {
            libMesh::Nemesis_IO io_l(*mesh_master);
            utopia::convert(lagr, *V_m.equation_system().solution);
            V_m.equation_system().solution->close();
            io_l.write_timestep("lagr.e", V_m.equation_systems(), 1, 0);
        } else {
            libMesh::Nemesis_IO io_l(*mesh_slave);
            utopia::convert(lagr, *V_s.equation_system().solution);
            V_s.equation_system().solution->close();
            io_l.write_timestep("lagr.e", V_s.equation_systems(), 1, 0);
        }
        
        c.stop();
        std::cout << c << std::endl;
    }
}

