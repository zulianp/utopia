#include "utopia_RMTRApp.hpp"

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
#include "utopia_MinSurf.hpp"
#include "utopia_Bratu.hpp"
#include "utopia_Poisson.hpp"
#include "utopia_IPTransfer.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"

namespace utopia {


    typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

    void RMTRApp::init(libMesh::Parallel::Communicator &comm)
    {
        comm_ = make_ref(comm);
    }

    static bool assemble_projection(FunctionSpaceT &from, FunctionSpaceT &to, USparseMatrix &T)
    {
        auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), true, false);
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

        UVector d = sum(*mats[0], 1);
        T = diag(1./d) * (*mats[0]);
        double sum_T = sum(T);
        std::cout << sum_T << std::endl;
        return true;
    }

    static bool assemble_interpolation(FunctionSpaceT &from, FunctionSpaceT &to, USparseMatrix &T)
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

        T = std::move(*mats[0]);
        return true;
    }


    static void refine(const int n_refs, libMesh::MeshBase &mesh)
    {
        if(n_refs <= 0) return;

        libMesh::MeshRefinement mesh_refinement(mesh);
        mesh_refinement.make_flags_parallel_consistent();
        mesh_refinement.uniformly_refine(n_refs);
    }

    class RMTRApp::SimulationInput : public Configurable {
    public:

        void read(Input &is) override
        {
            try {
                is.get("mesh", mesh_type);
                path = "";
                is.get("path", path);

                is.get("boundary-conditions", [this](Input &is) {
                    is.get_all([this](Input &is) {
                        int side_set = 0;

                        is.get("side", side_set);


                        double value = 0;
                        is.get("value", value);

                        sides.push_back(side_set);
                        values.push_back(value);
                    });
                });

                refinements = 0;
                is.get("refinements", refinements);

                order = 1;
                is.get("order", order);


                elem_type = "quad";
                is.get("elem-type", elem_type);

                fun = "bratu";
                is.get("function", fun);

                n_levels = 2;
                is.get("n-levels", n_levels);

                verbose = false;
                is.get("verbose", verbose);

                use_newton = false;
                is.get("use-newton", use_newton);

            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
            }
        }

        libMesh::ElemType get_type(const int dim) const
        {
            if(dim == 3) {
                libMesh::ElemType type = libMesh::HEX8;

                if(elem_type == "tet") {
                    type = libMesh::TET4;
                }

                if(order == 2) {
                    type = libMesh::HEX20;

                    if(elem_type == "tet") {
                        type = libMesh::TET10;
                    }
                }

                return type;

            } else if(dim == 2) {
                libMesh::ElemType type = libMesh::QUAD4;

                if(elem_type == "tri") {
                    type = libMesh::TRI3;
                }

                if(order == 2) {
                    type = libMesh::QUAD8;

                    if(elem_type == "tri") {
                        type = libMesh::TRI6;
                    }
                }

                return type;
            }

            return libMesh::TRI3;
        }

        void make_mesh(libMesh::DistributedMesh &mesh) const
        {
            if(this->mesh_type == "file") {
                mesh.read(path);
            } else if(this->mesh_type == "unit-square") {
                libMesh::MeshTools::Generation::build_square(mesh,
                                                             3, 3,
                                                             -0., 1.,
                                                             -0., 1.,
                                                             get_type(2)
                                                             );
            } else if(this->mesh_type == "aabb") {
                libMesh::DistributedMesh temp_mesh(mesh.comm());
                temp_mesh.read(path);

#if LIBMESH_VERSION_LESS_THAN(1, 3, 0)
                libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::bounding_box(temp_mesh);
#else
                libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::create_bounding_box(temp_mesh);
#endif

                if(temp_mesh.spatial_dimension() == 3) {

                    libMesh::MeshTools::Generation::build_cube(
                        mesh,
                        n, n, n,
                        bb.min()(0) - span[0], bb.max()(0) + span[0],
                        bb.min()(1) - span[1], bb.max()(1) + span[1],
                        bb.min()(2) - span[2], bb.max()(2) + span[2],
                        get_type(3)
                    );

                } else {

                    libMesh::MeshTools::Generation::build_square(
                        mesh,
                        n, n,
                        bb.min()(0) - span[0], bb.max()(0) + span[0],
                        bb.min()(1) - span[1], bb.max()(1) + span[1],
                        get_type(2)
                    );
                }

            }

            refine(refinements, mesh);

            if(this->mesh_type == "file" && order == 2) {
                mesh.all_second_order();
            }
        }

        void set_up_bc(FunctionSpaceT &V) const
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

        inline bool empty() const
        {
            return mesh_type.empty();
        }

        void describe(std::ostream &os = std::cout) const
        {
            os << "-----------------------------------\n";
            os << "mesh_type:       " << mesh_type << "\n";
            os << "order:           " << order << "\n";
            os << "path:            " << path << "\n";
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
        int refinements;
        int order;
        double span[3];
        std::string elem_type;
        int n = 10;
        std::string fun;
        int n_levels;
        bool verbose;
        bool use_newton;
    };

    static std::shared_ptr<ExtendedFunction<USparseMatrix, UVector>> get_function(const RMTRApp::SimulationInput &in, FunctionSpaceT &V)
    {
        std::shared_ptr<ExtendedFunction<USparseMatrix, UVector>> f;
        if(in.fun == "bratu") {
            f = std::make_shared<Bratu<decltype(V), USparseMatrix, UVector>>(V);
        } else if(in.fun == "min-surf") {
            f = std::make_shared<MinSurf<decltype(V), USparseMatrix, UVector>>(V);
        } else if(in.fun == "matrixpoisson") {
            f = std::make_shared<Poisson<decltype(V), USparseMatrix, UVector>>(V);
        } else if(in.fun == "poisson") {
            f = std::make_shared<FormPoisson<decltype(V), USparseMatrix, UVector>>(V);
        } 

        else {
            assert(false);
            return nullptr;
        }

        return f;
    }

    void RMTRApp::solve_newton(const SimulationInput &in)
    {
        Chrono c;
        c.start();

        auto mesh = std::make_shared<libMesh::DistributedMesh>(*comm_);
        in.make_mesh(*mesh);

        mesh->print_info();

        const auto elem_order = libMesh::Order(in.order);
        auto V = FunctionSpaceT(*mesh, libMesh::LAGRANGE, elem_order, "u");
        in.set_up_bc(V);
        V.initialize();
        std::cout << "n_dofs: " << V.dof_map().n_dofs() << std::endl;

        auto f = get_function(in, V);
        // auto f = std::make_shared<Poisson<decltype(V), USparseMatrix, UVector>>(V);

        Newton<USparseMatrix, UVector> solver;
        solver.set_linear_solver(std::make_shared<Factorization<USparseMatrix, UVector>>());
        // solver.set_line_search_strategy(std::make_shared<Backtracking<USparseMatrix, UVector>>());

        UVector x = local_zeros(V.dof_map().n_local_dofs());

        apply_boundary_conditions(V.dof_map(), x);
        solver.verbose(in.verbose);

        solver.solve(*f, x);

        double energy = 0.;
        f->value(x, energy);
        disp(energy);
        c.stop();
        std::cout << c << std::endl;

        write("rmtr.e", V, x);
    }


    void RMTRApp::solve_rmtr(const SimulationInput &in)
    {
        using TransferT   = utopia::Transfer<USparseMatrix, UVector>;
        using IPTransferT = utopia::IPTransfer<USparseMatrix, UVector>;

        std::size_t n_levels = in.n_levels;
        std::vector< std::shared_ptr<libMesh::DistributedMesh> > meshes(n_levels);
        std::vector< std::shared_ptr<FunctionSpaceT> > spaces(n_levels);
        std::vector< std::shared_ptr<ExtendedFunction<USparseMatrix, UVector>> > functions(n_levels);

        auto coarse_solver = std::make_shared<utopia::SteihaugToint<USparseMatrix, UVector, HOMEMADE> >();
        auto smoother      = std::make_shared<utopia::SteihaugToint<USparseMatrix, UVector, HOMEMADE> >();

        // coarse_solver->verbose(true); 
        // smoother->verbose(true); 
        // auto coarse_solver = std::make_shared<utopia::KSP_TR<DSMatrixd, DVectord> >("gltr");
        // coarse_solver->atol(1e-12);
        // coarse_solver->rtol(1e-12);
        // coarse_solver->pc_type("lu");



        meshes[0] = std::make_shared<libMesh::DistributedMesh>(*comm_);
        in.make_mesh(*meshes[0]);

        for(std::size_t i = 1; i < n_levels; ++i) {
            meshes[i] = std::make_shared<libMesh::DistributedMesh>(*meshes[i-1]);
            refine(1, *meshes[i]);
        }

        const auto elem_order = libMesh::Order(in.order);
        for(std::size_t i = 0; i < n_levels; ++i) {
            spaces[i]           = std::make_shared<FunctionSpaceT>(*meshes[i], libMesh::LAGRANGE, elem_order, "u");

            in.set_up_bc(*spaces[i]);
            spaces[i]->initialize();
        }

        std::vector< std::shared_ptr<TransferT> > transfers(n_levels - 1);

        for(std::size_t i = 1; i < n_levels; ++i) {
            auto T_cf = std::make_shared<USparseMatrix>();
            auto T_fc = std::make_shared<USparseMatrix>();
            assemble_interpolation(*spaces[i-1], *spaces[i],   *T_cf); //assemble_projection 
            assemble_interpolation(*spaces[i],   *spaces[i-1], *T_fc); //assemble_projection 
            transfers[i-1] = std::make_shared<IPTransferT>(T_cf, T_fc);
        }

        for(std::size_t i = 0; i < n_levels; ++i) {
            functions[i] = get_function(in, *spaces[i]);
        }

        auto rmtr = std::make_shared<RMTR<USparseMatrix, UVector, FIRST_ORDER> >(n_levels);
        rmtr->set_transfer_operators(transfers);

        
        rmtr->set_coarse_tr_strategy(coarse_solver); 
        rmtr->set_fine_tr_strategy(smoother); 

        rmtr->max_it(30);
        rmtr->max_coarse_it(3);
        rmtr->max_smoothing_it(3);
        rmtr->delta0(1000);
        rmtr->atol(1e-6);
        rmtr->rtol(1e-10);
        rmtr->set_grad_smoothess_termination(0.000001);
        rmtr->set_eps_grad_termination(1e-7);

        rmtr->verbose(in.verbose);
        // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
        rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
        rmtr->set_functions(functions);

        auto &dof_map = spaces.back()->dof_map();

        UVector x;
        rmtr->handle_equality_constraints();

        bool ok = rmtr->solve(x);

        std::cout<<"fine dofs:  "<< size(x).get(0) << "  \n"; 

        //Write solution to disk
        write("rmtr.e", *spaces.back(), x);
    }


    void RMTRApp::run(Input &in)
    {

        SimulationInput sim_in;
        in.get("rmtr-app", sim_in);

        sim_in.describe();

        if(sim_in.use_newton) {
            solve_newton(sim_in);
        } else {
            solve_rmtr(sim_in);
        }
    }
}

