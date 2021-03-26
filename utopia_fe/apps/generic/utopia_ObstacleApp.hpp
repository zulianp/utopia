#ifndef UTOPIA_OBSTACLE_APP
#define UTOPIA_OBSTACLE_APP

#include "utopia_fe_base.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_ILU.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class ObstacleProblem : public Configurable {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;

        // Use specialized compoenents for function space
        using Obstacle_t = utopia::Obstacle<FunctionSpace>;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;

        // Use algorithms from utopia algebra
        using QPSolver_t = utopia::QPSolver<Matrix_t, Vector_t>;
        using SemismoothNewton_t = utopia::SemismoothNewton<Matrix_t, Vector_t>;
        using Factorization_t = utopia::Factorization<Matrix_t, Vector_t>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using IterativeSolver_t = utopia::IterativeSolver<Matrix_t, Vector_t>;
        using AlgebraicMultigrid_t = utopia::AlgebraicMultigrid<Matrix_t, Vector_t>;
        using ProjectedGaussSeidel_t = utopia::ProjectedGaussSeidel<Matrix_t, Vector_t>;
        using KSPSolver_t = utopia::KSPSolver<Matrix_t, Vector_t>;

        void read(Input &in) override {
            in.get("space", [this](Input &in) {
                bool read_state = false;
                in.get("read_state", read_state);

                if (read_state) {
                    space.read_with_state(in, deformation);
                } else {
                    space.read(in);
                }
            });

            if (space.empty()) {
                valid = false;
                return;
            }

            in.get("obstacle", obstacle_mesh);
            in.get("obstacle", params);

            assembler = std::make_shared<OmniAssembler_t>(make_ref(space));
            in.get("assembly", *assembler);

            int block_size = space.mesh().spatial_dimension();
            bool use_amg = true;

            in.get("obstacle", [this](Input &in) {
                in.get("nonlinear_obstacle_iterations", nonlinear_contact_iterations);
                in.get("update_factor", update_factor);
            });

            std::string qp_solver_type = "pgs";
            in.get("qp_solver", [&qp_solver_type, &use_amg](Input &in) {
                in.get("use_amg", use_amg);
                in.get("type", qp_solver_type);
            });

            if (qp_solver_type == "ssnewton") {
                if (use_amg) {
                    std::shared_ptr<IterativeSolver_t> smoother;
                    auto ksp = std::make_shared<KSPSolver_t>();
                    ksp->pc_type("ilu");
                    ksp->ksp_type("richardson");
                    smoother = ksp;

                    std::shared_ptr<LinearSolver_t> coarse_solver;
                    coarse_solver = std::make_shared<KSPSolver_t>();

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

                    auto amg = std::make_shared<AlgebraicMultigrid_t>(smoother, coarse_solver, agglomerator);

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
                    auto linear_solver = std::make_shared<KSPSolver_t>();
                    qp_solver_ = std::make_shared<SemismoothNewton_t>(linear_solver);
                }

            } else {
                qp_solver_ = std::make_shared<ProjectedGaussSeidel_t>();
            }

            in.get("qp_solver", *qp_solver_);
        }

        void run() {
            obs.set_params(params);
            obs.init_obstacle(obstacle_mesh);
            obstacle_mesh.write("obstacle.e");

            Vector_t x, g, increment;
            Matrix_t H;

            space.create_matrix(H);
            space.create_vector(x);
            space.create_vector(g);

            if (!empty(deformation)) {
                x = deformation;
                space.displace(x);
            } else {
                x.set(0.0);
            }

            for (int i = 0; i < nonlinear_contact_iterations; ++i) {
                if (!obs.assemble(space)) {
                    Utopia::Abort();
                    return;
                }

                space.displace(-x);

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
                obs.inverse_transform(x_c, increment);
                x += update_factor * increment;

                // Output solution
                { space.write("obs_sol_" + std::to_string(i) + ".e", x); }

                space.displace(x);

                // Output extras
                {
                    Vector_t gap_zeroed = e_mul(obs.is_contact(), obs.gap());

                    // space.write("obs_gap.e", gap_zeroed);
                    // space.write("obs_is_contact.e", obs.is_contact());
                    // space.write("obs_normals.e", obs.normals());
                    // space.write("rhs.e", g);

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

                    space.write("rays_" + std::to_string(i) + ".e", rays);
                }
            }
        }

        ObstacleProblem() {}

        bool valid{true};

    private:
        FunctionSpace space;
        Vector_t deformation;

        Mesh_t obstacle_mesh;
        typename Obstacle_t::Params params;
        Obstacle_t obs;
        std::shared_ptr<OmniAssembler_t> assembler;
        std::shared_ptr<QPSolver_t> qp_solver_;

        int nonlinear_contact_iterations{1};
        Scalar_t update_factor{1};
    };

}  // namespace utopia

#endif
