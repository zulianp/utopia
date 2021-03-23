
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
#include "utopia_libmesh_OmniAssembler.hpp"

namespace utopia {

    namespace frontend = utopia::libmesh;

    class ObstacleProblem : public Configurable {
    public:
        using Vector_t = Traits<frontend::FunctionSpace>::Vector;
        using Matrix_t = Traits<frontend::FunctionSpace>::Matrix;
        using Size_t = Traits<frontend::FunctionSpace>::SizeType;
        using QPSolver_t = utopia::QPSolver<Matrix_t, Vector_t>;
        using SemismoothNewton_t = utopia::SemismoothNewton<Matrix_t, Vector_t>;
        using Factorization_t = utopia::Factorization<Matrix_t, Vector_t>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;

        void read(Input &in) override {
            in.get("space", space);

            if (space.empty()) {
                valid = false;
                return;
            }

            in.get("obstacle", obstacle_mesh);
            in.get("obstacle", params);

            assembler = std::make_shared<libmesh::OmniAssembler>(make_ref(space));
            in.get("assembly", *assembler);

            int block_size = space.mesh().spatial_dimension();
            bool use_amg = true;

            std::string qp_solver_type = "pgs";
            in.get("qp_solver", [&qp_solver_type, &use_amg](Input &in) {
                in.get("use_amg", use_amg);
                in.get("type", qp_solver_type);
            });

            if (qp_solver_type == "ssnewton") {
                if (use_amg) {
                    std::shared_ptr<IterativeSolver<Matrix_t, Vector_t>> smoother;
                    auto ksp = std::make_shared<KSPSolver<Matrix_t, Vector_t>>();
                    ksp->pc_type("ilu");
                    ksp->ksp_type("richardson");
                    smoother = ksp;

                    std::shared_ptr<LinearSolver_t> coarse_solver;
                    coarse_solver = std::make_shared<Factorization_t>();

                    std::shared_ptr<MatrixAgglomerator<Matrix_t>> agglomerator;

                    if (block_size == 2) {
                        agglomerator = std::make_shared<BlockAgglomerate<Matrix_t, 2>>();
                    } else if (block_size == 3) {
                        agglomerator = std::make_shared<BlockAgglomerate<Matrix_t, 3>>();
                    } else if (block_size == 4) {
                        agglomerator = std::make_shared<BlockAgglomerate<Matrix_t, 4>>();
                    } else {
                        agglomerator = std::make_shared<Agglomerate<Matrix_t>>();
                    }

                    auto amg =
                        std::make_shared<AlgebraicMultigrid<Matrix_t, Vector_t>>(smoother, coarse_solver, agglomerator);

                    InputParameters inner_params;
                    inner_params.set("block_size", block_size);
                    inner_params.set("use_simd", true);
                    inner_params.set("use_line_search", true);
                    coarse_solver->read(inner_params);
                    smoother->read(inner_params);
                    amg->verbose(true);

                    in.get("amg", *amg);

                    qp_solver_ = std::make_shared<SemismoothNewton_t>(amg);
                } else {
                    auto linear_solver = std::make_shared<Factorization_t>();
                    qp_solver_ = std::make_shared<SemismoothNewton_t>(linear_solver);
                }

            } else {
                qp_solver_ = std::make_shared<ProjectedGaussSeidel<Matrix_t, Vector_t>>();
            }

            in.get("qp_solver", *qp_solver_);
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

            Vector_t x, g;
            Matrix_t H;

            space.create_matrix(H);
            space.create_vector(x);
            space.create_vector(g);

            assembler->assemble(x, H, g);
            g *= -1.0;

            space.apply_constraints(H, g);

            Matrix_t H_c;
            Vector_t g_c(layout(x), 0.0), x_c(layout(x), 0.0);

            obs.transform(H, H_c);
            obs.transform(g, g_c);

            auto constr = make_upper_bound_constraints(make_ref(const_cast<Vector_t &>(obs.gap())));
            constr.fill_empty_bounds(layout(g_c));
            qp_solver_->set_box_constraints(constr);
            qp_solver_->solve(H_c, g_c, x_c);
            obs.inverse_transform(x_c, x);

            // Output solution
            { space.write("obs_sol.e", x); }

            // Output extras
            {
                obstacle_mesh.write("obstacle.e");
                Vector_t gap_zeroed = e_mul(obs.is_contact(), obs.gap());
                space.write("obs_gap.e", gap_zeroed);
                space.write("obs_is_contact.e", obs.is_contact());
                space.write("obs_normals.e", obs.normals());
                space.write("rhs.e", g);

                Vector_t rays = obs.normals();

                {
                    auto rays_view = local_view_device(rays);
                    auto gap_view = const_local_view_device(obs.gap());

                    int dim = space.mesh().spatial_dimension();

                    parallel_for(
                        local_range_device(rays), UTOPIA_LAMBDA(const Size_t i) {
                            Size_t block = i / dim;
                            Size_t g_idx = dim * block;

                            auto ri = rays_view.get(i);
                            ri *= gap_view.get(g_idx);
                            rays_view.set(i, ri);
                        });
                }

                space.write("rays.e", rays);
            }
        }

        ObstacleProblem() {}

        bool valid{true};

    private:
        frontend::FunctionSpace space;
        frontend::Mesh obstacle_mesh;
        frontend::Obstacle::Params params;
        frontend::Obstacle obs;
        std::shared_ptr<frontend::OmniAssembler> assembler;
        std::shared_ptr<QPSolver_t> qp_solver_;
    };

}  // namespace utopia

void obs(utopia::Input &in) {
    utopia::ObstacleProblem obs;
    obs.read(in);
    if (obs.valid) {
        obs.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(obs);
