#include "utopia_ObstacleApp.hpp"

#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
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

        void read(Input &in) override {
            in.get("mesh", mesh);
            in.get("space", space);

            in.get("obstacle", obstacle_mesh);
            in.get("obstacle", params);
        }

        void run() {
            obs.init(obstacle_mesh.mesh());
            obs.assemble(mesh.mesh(), space.space()[0].dof_map(), params);

            write("obstacle.e", obstacle_mesh.mesh());
            write("obs_sol.e", space.space()[0], obs.output().gap);
        }

        ObstacleProblem(libMesh::Parallel::Communicator &comm)
            : mesh(comm), space(make_ref(mesh)), obstacle_mesh(comm) {}

    private:
        UIMesh<libMesh::DistributedMesh> mesh;
        UIFunctionSpace<FS> space;

        UIMesh<libMesh::DistributedMesh> obstacle_mesh;

        ObstacleParams params;
        Obstacle obs;
    };

    void ObstacleApp::run(Input &in) {
        ObstacleProblem obs(comm());
        obs.read(in);
        obs.run();
    }
}  // namespace utopia
