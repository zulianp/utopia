#include "utopia_ObstacleApp.hpp"

#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMaterial.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_ui.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_ILU.hpp"

#include "utopia_MPITimeStatistics.hpp"
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_Obstacle.hpp"

#include "libmesh/boundary_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/nemesis_io.h"

namespace utopia {

    inline void write(const Path &path, libMesh::MeshBase &mesh) {
        std::cout << "Writing mesh at " << path.to_string() << std::endl;
        libMesh::ExodusII_IO(mesh).write(path.to_string());
    }

    class ObstacleProblem : public Configurable {
    public:
        using FS = utopia::LibMeshFunctionSpace;
        using PFS = utopia::ProductFunctionSpace<FS>;
        using MaterialT = utopia::UIMaterial<PFS, USparseMatrix, UVector>;
        using ForcingFunctionT = utopia::UIForcingFunction<PFS, UVector>;

        void read(Input &in) override {
            in.get("mesh", mesh);
            in.get("space", space);

            in.get("obstacle", obstacle_mesh);
            in.get("obstacle", params);

            auto model = make_unique<MaterialT>(space.space());
            auto forcing_function = make_unique<ForcingFunctionT>(space.space());

            in.get("model", *model);
            in.get("forcing-functions", *forcing_function);

            assert(model->good());

            model_ =
                std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(std::move(model), std::move(forcing_function));

            int block_size = mesh.mesh().spatial_dimension();
            bool use_amg = true;
            in.get("use_amg", use_amg);

            std::string qp_solver_type = "pgs";
            in.get("qp-solver", [&qp_solver_type](Input &in) { in.get("type", qp_solver_type); });

            if (qp_solver_type == "ssnewton") {
                if (use_amg) {
                    std::shared_ptr<IterativeSolver<USparseMatrix, UVector>> smoother;
                    auto ksp = std::make_shared<KSPSolver<USparseMatrix, UVector>>();
                    ksp->pc_type("ilu");
                    ksp->ksp_type("richardson");
                    smoother = ksp;

                    std::shared_ptr<LinearSolver<USparseMatrix, UVector>> coarse_solver;
                    coarse_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();

                    std::shared_ptr<MatrixAgglomerator<USparseMatrix>> agglomerator;

                    if (block_size == 2) {
                        agglomerator = std::make_shared<BlockAgglomerate<USparseMatrix, 2>>();
                    } else if (block_size == 3) {
                        agglomerator = std::make_shared<BlockAgglomerate<USparseMatrix, 3>>();
                    } else if (block_size == 4) {
                        agglomerator = std::make_shared<BlockAgglomerate<USparseMatrix, 4>>();
                    } else {
                        agglomerator = std::make_shared<Agglomerate<USparseMatrix>>();
                    }

                    auto amg = std::make_shared<AlgebraicMultigrid<USparseMatrix, UVector>>(
                        smoother, coarse_solver, agglomerator);

                    InputParameters inner_params;
                    inner_params.set("block_size", block_size);
                    inner_params.set("use_simd", true);
                    inner_params.set("use_line_search", true);
                    coarse_solver->read(inner_params);
                    smoother->read(inner_params);
                    amg->verbose(true);

                    in.get("amg", *amg);

                    qp_solver_ = std::make_shared<SemismoothNewton<USparseMatrix, UVector>>(amg);
                } else {
                    auto linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
                    qp_solver_ = std::make_shared<SemismoothNewton<USparseMatrix, UVector>>(linear_solver);
                }

            } else {
                qp_solver_ = std::make_shared<ProjectedGaussSeidel<USparseMatrix, UVector>>();
            }

            in.get("qp-solver", *qp_solver_);
        }

        void run() {
            obs.init(obstacle_mesh.mesh());
            obs.assemble(mesh.mesh(), space.space()[0].dof_map(), params);

            auto &dof_map = space.space()[0].dof_map();

            UVector x, g, c;
            USparseMatrix H;

            UIndexSet ghost_nodes;
            convert(dof_map.get_send_list(), ghost_nodes);
            x = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
            c = x;

            model_->assemble_hessian_and_gradient(x, H, g);
            g *= -1.0;
            apply_boundary_conditions(dof_map, H, g);

            USparseMatrix H_c;
            UVector g_c(layout(x), 0.0), x_c(layout(x), 0.0);

            obs.transform(H, H_c);
            obs.transform(g, g_c);

            auto constr = make_upper_bound_constraints(make_ref(obs.output().gap));
            constr.fill_empty_bounds(layout(g_c));
            qp_solver_->set_box_constraints(constr);

            // Output
            {
                write("obstacle.e", obstacle_mesh.mesh());
                UVector gap_zeroed = e_mul(obs.output().is_contact, obs.output().gap);
                write("obs_gap.e", space.space()[0], gap_zeroed);
                write("obs_normals.e", space.space()[0], obs.output().normals);
                write("obs_is_contact.e", space.space()[0], obs.output().is_contact);
            }

            qp_solver_->solve(H_c, g_c, x_c);

            obs.inverse_transform(x_c, x);

            // Output
            { write("obs_sol.e", space.space()[0], x); }
        }

        ObstacleProblem(libMesh::Parallel::Communicator &comm)
            : mesh(comm), space(make_ref(mesh)), obstacle_mesh(comm) {}

    private:
        UIMesh<libMesh::DistributedMesh> mesh;
        UIFunctionSpace<FS> space;

        UIMesh<libMesh::DistributedMesh> obstacle_mesh;

        ObstacleParams params;
        Obstacle obs;

        std::shared_ptr<ElasticMaterial<USparseMatrix, UVector>> model_;
        std::shared_ptr<QPSolver<USparseMatrix, UVector>> qp_solver_;
    };

    void ObstacleApp::run(Input &in) {
        ObstacleProblem obs(comm());
        obs.read(in);
        obs.run();
    }
}  // namespace utopia
