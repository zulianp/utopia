
#include "utopia_Main.hpp"

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

#include "utopia_libmesh_Obstacle.hpp"

#include "libmesh/boundary_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/nemesis_io.h"

namespace utopia {

    class ObstacleProblem : public Configurable {
    public:
        using FS = utopia::LibMeshFunctionSpace;
        using PFS = utopia::ProductFunctionSpace<FS>;
        using MaterialT = utopia::UIMaterial<PFS, USparseMatrix, UVector>;
        using ForcingFunctionT = utopia::UIForcingFunction<PFS, UVector>;

        void read(Input &in) override {
            in.get("space", space);
            in.get("obstacle", obstacle_mesh);
            in.get("obstacle", params);

            // hack_ = std::make_shared<LibMeshFunctionSpace>(space.raw_type(), 0);
            // auto model = make_unique<MaterialT>(*hack_);
            // auto forcing_function = make_unique<ForcingFunctionT>(space.space());
            // in.get("model", *model);
            // in.get("forcing-functions", *forcing_function);

            // assert(model->good());

            // model_ =
            //     std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(std::move(model),
            //     std::move(forcing_function));

            // model_ = std::move(model);

            int block_size = space.mesh().spatial_dimension();
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
            obs.set_params(params);
            obs.init_obstacle(obstacle_mesh);

            if (!obs.assemble(space)) {
                Utopia::Abort();
                return;
            }

            // FIXME
            auto &dof_map = space.raw_type_dof_map();

            UVector x, g;
            USparseMatrix H;

            space.create_matrix(H);
            space.create_vector(x);
            space.create_vector(g);

            model_->assemble_hessian_and_gradient(x, H, g);
            g *= -1.0;

            space.apply_constraints(H, g);

            USparseMatrix H_c;
            UVector g_c(layout(x), 0.0), x_c(layout(x), 0.0);

            obs.transform(H, H_c);
            obs.transform(g, g_c);

            auto constr = make_upper_bound_constraints(make_ref(const_cast<UVector &>(obs.gap())));
            constr.fill_empty_bounds(layout(g_c));
            qp_solver_->set_box_constraints(constr);

            // Output
            {
                obstacle_mesh.write("obstacle.e");
                UVector gap_zeroed = e_mul(obs.is_contact(), obs.gap());
                space.write("obs_gap.e", gap_zeroed);
                space.write("obs_is_contact.e", obs.is_contact());
            }

            qp_solver_->solve(H_c, g_c, x_c);

            obs.inverse_transform(x_c, x);

            // Output
            { space.write("obs_sol.e", x); }
        }

        ObstacleProblem() {}

    private:
        utopia::libmesh::FunctionSpace space;
        utopia::libmesh::Mesh obstacle_mesh;
        utopia::libmesh::Obstacle::Params params;
        utopia::libmesh::Obstacle obs;

        std::shared_ptr<ElasticMaterial<USparseMatrix, UVector>> model_;
        std::shared_ptr<QPSolver<USparseMatrix, UVector>> qp_solver_;
    };

}  // namespace utopia

void obs(utopia::Input &in) {
    utopia::ObstacleProblem obs;
    obs.read(in);
    obs.run();
}

UTOPIA_REGISTER_APP(obs);
