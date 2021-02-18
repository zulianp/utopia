#include "utopia_ObstacleApp.hpp"

#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMaterial.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_ui.hpp"

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

            qp_solver_ = std::make_shared<ProjectedGaussSeidel<USparseMatrix, UVector>>();
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

            qp_solver_->solve(H_c, g_c, x_c);

            obs.inverse_transform(x_c, x);

            write("obstacle.e", obstacle_mesh.mesh());
            write("obs_gap.e", space.space()[0], obs.output().gap);
            write("obs_normals.e", space.space()[0], obs.output().normals);
            write("obs_is_contact.e", space.space()[0], obs.output().is_contact);
            write("obs_sol.e", space.space()[0], x);
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
