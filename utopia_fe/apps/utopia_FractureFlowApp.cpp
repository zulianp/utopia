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



#include "libmesh/mesh_refinement.h"

namespace utopia {
    
    
    typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
    
    void FractureFlowApp::init(libMesh::LibMeshInit &init)
    {
        comm_ = make_ref(init.comm());
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
        return true;
    }
    
    
    static void solve_monolithic(
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
        
        DSMatrixd B, D;
        assemble_projection(V_m, V_s, B, D);
        
        D *= -1.;
        // B *= -1.;
        
        auto s_m = local_size(A_m);
        auto s_s = local_size(A_s);
        
        
        auto nnz_x_row_m =
        std::max(*std::max_element(V_m.dof_map().get_n_nz().begin(), V_m.dof_map().get_n_nz().end()),
                 *std::max_element(V_m.dof_map().get_n_oz().begin(), V_m.dof_map().get_n_oz().end()));
        
        auto nnz_x_row_s =
        std::max(*std::max_element(V_s.dof_map().get_n_nz().begin(), V_s.dof_map().get_n_nz().end()),
                 *std::max_element(V_s.dof_map().get_n_oz().begin(), V_s.dof_map().get_n_oz().end()));
        
        
        DSMatrixd A = local_sparse(
                                   s_m.get(0) + 2 * s_s.get(0),
                                   s_m.get(1) + 2 * s_s.get(1),
                                   2 * nnz_x_row_m + 2 * nnz_x_row_s //FIXME
                                   );
        
        
        auto rr_m = row_range(A_m);
        auto rr_s = row_range(A_s);
        
        auto off_r   = size(A_m).get(0);
        auto off_c   = size(A_m).get(1);
        
        auto off_r_l = off_r + size(A_s).get(0);
        auto off_c_l = off_c + size(A_s).get(1);
        
        {
            Write<DSMatrixd> w_(A);
            
            for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
                {
                    RowView<DSMatrixd> row(A_m, i);
                    
                    for(auto k = 0; k < row.n_values(); ++k) {
                        //(1, 1)
                        A.set(i, row.col(k), row.get(k));
                    }
                }
            }
            
            for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
                
                {
                    RowView<DSMatrixd> row(A_s, i);
                    
                    for(auto k = 0; k < row.n_values(); ++k) {
                        //(2, 2)
                        A.set(off_r + i, off_c + row.col(k), row.get(k));
                    }
                }
                
                {
                    RowView<DSMatrixd> row(B, i);
                    
                    for(auto k = 0; k < row.n_values(); ++k) {
                        //(1, 3)  //B^T
                        if(!V_m.dof_map().is_constrained_dof(row.col(k))) {
                            A.set(row.col(k), off_c_l + i, row.get(k));
                        }
                        
                        //(3, 1) //B
                        A.set(off_r_l + i, row.col(k), row.get(k));
                    }
                }
                
                {
                    RowView<DSMatrixd> row(D, i);
                    
                    for(auto k = 0; k < row.n_values(); ++k) {
                        
                        //(2, 3) //D^T
                        if(!V_s.dof_map().is_constrained_dof(row.col(k))) {
                            A.set(off_r + row.col(k), off_c_l + i, row.get(k));
                        }
                        
                        //(3, 2)   //D
                        A.set(off_r_l + i, off_c + row.col(k), row.get(k));
                    }
                }
                
            }
        }
        
        DVectord rhs = local_zeros(local_size(rhs_m).get(0) + 2 * local_size(rhs_s).get(0));
        
        {
            Write<DVectord> w_(rhs);
            Write<DVectord> r_m(rhs_m), r_s(rhs_s);
            
            for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
                rhs.set(i, rhs_m.get(i));
            }
            
            for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
                rhs.set(i, off_r + rhs_s.get(i));
            }
        }
        
        // A.implementation().set_name("a");
        // write("A.m", A);
        
        // Factorization<DSMatrixd, DVectord> op;
        LUDecomposition<DSMatrixd, DVectord> op;
        op.update(make_ref(A));
        
        
        DVectord sol = local_zeros(local_size(rhs));
        op.apply(rhs, sol);
        
        
        sol_m  = local_zeros(local_size(rhs_m));
        sol_s  = local_zeros(local_size(rhs_s));
        lagr   = local_zeros(local_size(rhs_s));
        
        {
            Write<DVectord> w_(sol);
            Write<DVectord> r_m(sol_m), r_s(sol_s);
            
            for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
                sol_m.set(i, sol.get(i));
            }
            
            for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
                sol_s.set(i, sol.get(i) - off_r);
            }
            
            for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
                lagr.set(i, sol.get(i) - off_r_l);
            }
        }
        
    }
    
    static void write_solution(
                               const std::string &name,
                               DVectord &sol,
                               const int time_step,
                               const double t,
                               FunctionSpaceT &space,
                               libMesh::ExodusII_IO &io)
    {
        utopia::convert(sol, *space.equation_system().solution);
        space.equation_system().solution->close();
        io.write_timestep(name, space.equation_systems(), time_step, t);
    }
    
    static void solve_staggered(
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
        
        libMesh::ExodusII_IO io_m(V_m.mesh());
        libMesh::ExodusII_IO io_s(V_s.mesh());
        libMesh::ExodusII_IO io_m_l(V_m.mesh());
        libMesh::ExodusII_IO io_s_l(V_s.mesh());
        
        
        Factorization<DSMatrixd, DVectord> op_m;
        op_m.update(make_ref(A_m));
        
        
        Factorization<DSMatrixd, DVectord> op_s;
        op_s.update(make_ref(A_s));
        
        sol_s = local_zeros(local_size(rhs_s));
        
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
        
        
        // t.initialize(INTERPOLATION);
        // t.write("T.m");
        t.initialize(L2_PROJECTION);
        
        lagr_m.set(0.);
        
        for(int i = 0; i < 20; ++i) {
            apply_zero_boundary_conditions(V_m.dof_map(), lagr_m);
            rhs_lagr_m = rhs_m + lagr_m;
            
            write_solution(
                           "lagr_m.e",
                           lagr_m,
                           (i + 1),
                           (i + 1),
                           V_m,
                           io_m_l
                           );
            
            op_m.apply(rhs_lagr_m, sol_m);
            
            write_solution(
                           "sol_m.e",
                           sol_m,
                           (i + 1),
                           (i + 1),
                           V_m,
                           io_m
                           );
            
            sol_s.set(0.);
            t.apply(sol_m, sol_s);
            
            // apply_boundary_conditions(V_s, A_s, sol_s);
            
            lagr_s = rhs_s - A_s * sol_s;
            apply_zero_boundary_conditions(V_s.dof_map(), lagr_s);
            
            // rhs_lagr_s = rhs_s + lagr_s;
            // op_s.apply(rhs_lagr_s, sol_s);
            
            // lagr_s = rhs_lagr_s - A_s * sol_s;
            
            write_solution(
                           "sol_s.e",
                           sol_s,
                           (i + 1),
                           (i + 1),
                           V_s,
                           io_s
                           );
            
            
            write_solution(
                           "lagr_s.e",
                           lagr_s,
                           (i + 1),
                           (i + 1),
                           V_s,
                           io_s_l
                           );
            
            apply_zero_boundary_conditions(V_s.dof_map(), lagr_s);
            
            double n_lagr_s = norm2(lagr_s);
            
            disp(n_lagr_s);
            
            if(n_lagr_s < 1e-14) {
                break;
            }
            
            delta_lagr.set(0);
            t.apply_transpose(lagr_s, delta_lagr);
            // lagr_m += dumping * delta_lagr;
            lagr_m =  delta_lagr;
        }
        
        lagr = lagr_m;
        
    }
    
    static void solve_separate(
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
        
        t.initialize(INTERPOLATION);
        
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
        }

        void make_mesh(libMesh::DistributedMesh &mesh) const
        {
            if(this->mesh_type == "file") {
                mesh.read(path);
            } else if(this->mesh_type == "unit-square") {
                const auto elem_type = libMesh::QUAD8;
                libMesh::MeshTools::Generation::build_square(mesh,
                                                             30, 30,
                                                             -0., 1.,
                                                             -0., 1.,
                                                             elem_type
                                                             );
            }
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

        void desribe(std::ostream &os = std::cout) const
        {
            os << "mesh_type: " << mesh_type << "\n";
            os << "path: " << path << "\n";

            os << "side, value\n";

            std::size_t n = sides.size();
            for(std::size_t i = 0; i < n; ++i) {
                os << sides[i]<< ", " << values[i] << "\n";
            }   
        }
        

        std::string mesh_type, path;
        std::vector<int> sides;
        std::vector<double> values;
    };
    
    void FractureFlowApp::run(const std::string &conf_file_path)
    {
        Chrono c;
        
        c.start();
        
        auto is_ptr = open_istream(conf_file_path);

        Input master_in, slave_in;

        is_ptr->read("master", master_in);
        is_ptr->read("slave", slave_in);
        
        std::string solve_strategy = "staggered";
        is_ptr->read("solve-strategy", solve_strategy);

        master_in.desribe();
        slave_in.desribe();

        std::cout << "solve_strategy: "  << solve_strategy << std::endl;
        //model parameters
        // const unsigned int n = 100;
        // const unsigned int m = 10;
        
        //discretization parameters
        // const auto elem_type = libMesh::QUAD8;
        const auto elem_order = libMesh::FIRST;
        
        auto mesh_master = std::make_shared<libMesh::DistributedMesh>(*comm_);
        auto mesh_slave = std::make_shared<libMesh::DistributedMesh>(*comm_);

        master_in.make_mesh(*mesh_master);
        slave_in.make_mesh(*mesh_slave);



        // if( )
        // libMesh::MeshTools::Generation::build_square(*mesh_master,
        //                                              30, 30,
        //                                              // -0.5, 0.5,
        //                                              // -0.1, 0.1,
        //                                              -0., 1.,
        //                                              -0., 1.,
        //                                              elem_type
        //                                              );
        
        
        // mesh_master->read("../data/frac/master_backg.e");
        
        
        // libMesh::MeshTools::Generation::build_square(
        // 	*mesh_slave,
        // 	m, m/2,
        // 	0., 1.,
        // 	0.35, 0.4,
        // 	elem_type
        // );
        
        // mesh_slave->read("../data/frac/slave_frac.e");
        // mesh_slave->read("../data/frac/line.e");
        // refine(4, *mesh_slave);
        
        
        //equations system
        auto equation_systems_master = std::make_shared<libMesh::EquationSystems>(*mesh_master);
        auto &sys_master = equation_systems_master->add_system<libMesh::LinearImplicitSystem>("master");
        
        auto equation_systems_slave = std::make_shared<libMesh::EquationSystems>(*mesh_slave);
        auto &sys_slave = equation_systems_slave->add_system<libMesh::LinearImplicitSystem>("slave");
        
        //scalar function space
        auto V_m = FunctionSpaceT(equation_systems_master, libMesh::LAGRANGE, elem_order, "u");
        auto V_s = FunctionSpaceT(equation_systems_slave, libMesh::LAGRANGE, elem_order, "u");
        
        
        auto u_m = trial(V_m);
        auto v_m = test(V_m);
        
        auto u_s = trial(V_s);
        auto v_s = test(V_s);

        master_in.set_up_bc(V_m);
        slave_in.set_up_bc(V_s);
        
        V_m.initialize();
        V_s.initialize();
        
        auto eq_m = 1. * inner(grad(u_m), grad(v_m)) * dX == inner(coeff(0.), v_m) * dX;
        auto eq_s = 10. * inner(grad(u_s), grad(v_s)) * dX == inner(coeff(0.), v_s) * dX;
        
        DSMatrixd A_m, A_s;
        DVectord rhs_m, rhs_s;
        utopia::assemble(eq_m, A_m, rhs_m);
        utopia::assemble(eq_s, A_s, rhs_s);
        
        apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
        apply_boundary_conditions(V_s.dof_map(), A_s, rhs_s);
        
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        DVectord sol_m, sol_s, lagr;
        

        if(solve_strategy == "staggered") {
        // solve_monolithic(
        solve_staggered(
                        // solve_separate(
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
        solve_monolithic(
        // solve_staggered(
                        // solve_separate(
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
    }
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        libMesh::ExodusII_IO io_m(*mesh_master);
        libMesh::ExodusII_IO io_s(*mesh_slave);
        
        
        utopia::convert(sol_m, *V_m.equation_system().solution);
        V_m.equation_system().solution->close();
        io_m.write_timestep(V_m.equation_system().name() + ".e", V_m.equation_systems(), 1, 0);
        
        
        utopia::convert(sol_s, *V_s.equation_system().solution);
        V_s.equation_system().solution->close();
        io_s.write_timestep(V_s.equation_system().name() + ".e", V_s.equation_systems(), 1, 0);
        
        if(size(lagr) == size(sol_m)) {
            libMesh::ExodusII_IO io_l(*mesh_master);
            utopia::convert(lagr, *V_m.equation_system().solution);
            V_m.equation_system().solution->close();
            io_l.write_timestep("lagr.e", V_m.equation_systems(), 1, 0);
        } else {
            libMesh::ExodusII_IO io_l(*mesh_slave);
            utopia::convert(lagr, *V_s.equation_system().solution);
            V_s.equation_system().solution->close();
            io_l.write_timestep("lagr.e", V_s.equation_systems(), 1, 0);
        }
        
        
        
        
        
        c.stop();
        std::cout << c << std::endl;
    }
}

