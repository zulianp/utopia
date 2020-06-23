#include "utopia_TransferApp.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_ui.hpp"

#include "libmesh/boundary_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/node.h"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"

namespace utopia {

    class TransferApp::InputSpace : public Configurable {
    public:
        InputSpace(libMesh::Parallel::Communicator &comm) : mesh_(comm), space_(make_ref(mesh_)) {}

        void read(Input &is) override {
            try {
                is.get("mesh", mesh_);
                is.get("space", space_);

            } catch (const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
            }
        }

        inline bool empty() const { return mesh_.empty(); }

        inline libMesh::MeshBase &mesh() { return mesh_.mesh(); }

        inline std::shared_ptr<libMesh::MeshBase> mesh_ptr() { return mesh_.mesh_ptr(); }

        inline LibMeshFunctionSpace &space() { return space_.space()[0]; }

        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<LibMeshFunctionSpace> space_;
    };

    void TransferApp::run(Input &in) {
        Chrono c;
        c.start();

        InputSpace input_master(comm());
        InputSpace input_slave(comm());
        bool master_boundary = false, slave_boundary = false;

        std::shared_ptr<libMesh::MeshBase> master_actual_mesh, slave_actual_mesh;
        std::shared_ptr<LibMeshFunctionSpace> master_actual_space, slave_actual_space;

        std::shared_ptr<MeshTransferOperator> transfer_operator;

        in.get("transfer", [&](Input &is) {
            // get spaces
            is.get("master", input_master);
            is.get("slave", input_slave);

            std::cout << "1) read fe spaces form disk" << std::endl;

            is.get("master-boundary", master_boundary);
            is.get("slave-boundary", slave_boundary);

            if (master_boundary) {
                auto b_mesh = std::make_shared<libMesh::BoundaryMesh>(comm(), input_master.mesh().mesh_dimension() - 1);
                input_master.mesh().boundary_info->sync(*b_mesh);
                master_actual_space =
                    std::make_shared<LibMeshFunctionSpace>(*b_mesh, libMesh::LAGRANGE, libMesh::FIRST, "u");
                master_actual_space->initialize();
                master_actual_mesh = b_mesh;
            } else {
                master_actual_mesh = input_master.mesh_ptr();
                master_actual_space = make_ref(input_master.space());
            }

            if (slave_boundary) {
                auto b_mesh = std::make_shared<libMesh::BoundaryMesh>(comm(), input_slave.mesh().mesh_dimension() - 1);
                input_slave.mesh().boundary_info->sync(*b_mesh);
                slave_actual_space =
                    std::make_shared<LibMeshFunctionSpace>(*b_mesh, libMesh::LAGRANGE, libMesh::FIRST, "u");
                slave_actual_space->initialize();
                slave_actual_mesh = b_mesh;
            } else {
                slave_actual_mesh = input_slave.mesh_ptr();
                slave_actual_space = make_ref(input_slave.space());
            }

            transfer_operator = std::make_shared<MeshTransferOperator>(master_actual_mesh,
                                                                       make_ref(master_actual_space->dof_map()),
                                                                       slave_actual_mesh,
                                                                       make_ref(slave_actual_space->dof_map()));

            transfer_operator->read(is);

#ifdef WITH_TINY_EXPR
            std::string expr = "x";
            is.get("function", expr);

            fun = std::make_shared<SymbolicFunction>(expr);

            std::string fun_type = "non-linear";
            is.get("function-type", fun_type);

            if (fun_type == "constant") {
                fun_is_constant = true;
            } else {
                fun_is_constant = false;
            }
#else
            double expr = 1.;
            is.get("function", expr);
            fun_is_constant = true;

            fun = std::make_shared<ConstantCoefficient<double, 0>>(expr);
#endif  // WITH_TINY_EXPR
        });

        c.stop();

        if (!transfer_operator) {
            std::cout << "[Error] unable to read setting file" << std::endl;
            return;
        }

        std::cout << "2) assembling transfer operator..." << std::endl;

        bool ok = transfer_operator->assemble();

        std::cout << "3) assembly of transfer operator DONE" << std::endl;

        if (!ok) {
            std::cerr << "[Error] unable to assemble operator" << std::endl;
            assert(false);
            exit(-1);
        }

        auto u = trial(input_master.space());
        auto v = test(input_master.space());

        c.start();

        UVector fun_master_h, fun_master, fun_slave, back_fun_master;

        USparseMatrix mass_mat_master;
        assemble(inner(u, v) * dX, mass_mat_master);

        if (!fun_is_constant) {
            assemble(inner(*fun, v) * dX, fun_master_h);

            fun_master = fun_master_h;

            BiCGStab<USparseMatrix, UVector> solver;
            solver.solve(mass_mat_master, fun_master_h, fun_master);
        } else {
#ifdef WITH_TINY_EXPR
            fun_master = local_values(input_master.space().dof_map().n_local_dofs(), fun->eval(0., 0., 0.));
#else
            fun_master = local_values(input_master.space().dof_map().n_local_dofs(), fun->expr());
#endif  // WITH_TINY_EXPR
        }

        c.stop();

        if (mpi_world_rank() == 0) {
            std::cout << "Assembled M and fun_m" << std::endl;
            std::cout << c << std::endl;
            std::cout << std::flush;
        }

        mpi_world_barrier();

        transfer_operator->apply(fun_master, fun_slave);
        transfer_operator->apply_transpose(fun_slave, back_fun_master);

        double sum_fun_master = sum(fun_master);
        double sum_fun_slave = sum(fun_slave);

        std::cout << "master: " << sum_fun_master << " slave: " << sum_fun_slave << std::endl;

        ////////////////////////////////////////////////////////////
        // output
        ////////////////////////////////////////////////////////////

        convert(fun_master, *input_master.space().equation_system().solution);
        input_master.space().equation_system().solution->close();

        libMesh::Nemesis_IO io_master(input_master.mesh());
        io_master.write_equation_systems("master.e", input_master.space().equation_systems());

        convert(back_fun_master, *input_master.space().equation_system().solution);
        input_master.space().equation_system().solution->close();

        libMesh::Nemesis_IO io_master_adj(input_master.mesh());
        io_master_adj.write_equation_systems("master_adj.e", input_master.space().equation_systems());

        ////////////////////////////////////////////////////////////
        convert(fun_slave, *input_slave.space().equation_system().solution);
        input_slave.space().equation_system().solution->close();

        if (mpi_world_size() == 1) {
            libMesh::ExodusII_IO io_slave(input_slave.mesh());
            io_slave.write_equation_systems("slave.e", input_slave.space().equation_systems());
        } else {
            libMesh::Nemesis_IO io_slave(input_slave.mesh());
            io_slave.write_equation_systems("slave.e", input_slave.space().equation_systems());
        }
    }

    TransferApp::TransferApp() {}
    TransferApp::~TransferApp() {}
}  // namespace utopia
