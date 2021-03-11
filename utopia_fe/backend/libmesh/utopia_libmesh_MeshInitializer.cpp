#include "utopia_libmesh_MeshInitializer.hpp"
#include "utopia_Options.hpp"

#include "utopia_fe_base.hpp"

// All libmesh includes
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h"

namespace utopia {
    namespace libmesh {

        class MeshInitializer::Build {
        public:
            static constexpr int max_dim = 3;

            void default_cube(unsigned int default_n = 10,
                              libMesh::Real default_min = 0.0,
                              libMesh::Real default_max = 1.0) {
                for (auto &ni : n) {
                    ni = default_n;
                }

                for (auto &min_i : min) {
                    min_i = default_min;
                }

                for (auto &max_i : max) {
                    max_i = default_max;
                }
            }

            void default_square(unsigned int default_n = 10,
                                libMesh::Real default_min = 0.0,
                                libMesh::Real default_max = 1.0) {
                n[0] = n[1] = default_n;
                n[2] = 0;

                min[0] = min[1] = default_min;
                max[0] = max[1] = default_max;
                min[2] = max[2] = 0.0;
            }

            void default_line(unsigned int default_n = 10,
                              libMesh::Real default_min = 0.0,
                              libMesh::Real default_max = 1.0) {
                n[0] = default_n;
                n[1] = n[2] = 0;

                min[0] = min[1] = default_min;
                min[1] = max[1] = 0.0;
                min[2] = max[2] = 0.0;
            }

            void zero() { default_cube(0, 0.0, 0.0); }

            int n[max_dim];
            libMesh::Real min[max_dim], max[max_dim];
        };

        MeshInitializer::MeshInitializer(Mesh &mesh) : mesh_(mesh) {}

        void MeshInitializer::read(Input &in) {
            std::string parallelism = "distributed";
            std::string type = "cube";
            std::string elem_type = "HEX8";
            int n = 10;
            libMesh::Real min = 0.0, max = 1.0;
            bool gausss_lobatto_grid = false;

            Build b;
            b.default_cube(n, min, max);

            if (!Options()
                     .add_option(
                         "parallelism", parallelism, "Parallelism avaialble in libMesh: distributed|replicated|serial.")
                     .add_option("elem_type", elem_type, "Mesh element type. See libMesh ElemType enum.")
                     .add_option("n", n, "Number of nodes for all dimensions.")
                     .add_option("min", min, "Minumum coordinate for all directions.")
                     .add_option("max", max, "Maximum coordinate for all directions.")
                     .add_option("gausss_lobatto_grid", gausss_lobatto_grid, "Use Gauss-Lobatto grid.")
                     .add_option("nx", b.n[0], "Number of nodes in x direction")
                     .add_option("ny", b.n[1], "Number of nodes in y direction")
                     .add_option("nz", b.n[2], "Number of nodes in z direction")
                     .add_option("min_x", b.n[0], "Minumum coordinate in x direction")
                     .add_option("min_y", b.n[1], "Minumum coordinate in y direction")
                     .add_option("min_z", b.n[2], "Minumum coordinate in z direction")
                     .add_option("max_x", b.n[0], "Maximum coordinate in x direction")
                     .add_option("max_y", b.n[1], "Maximum coordinate in y direction")
                     .add_option("max_z", b.n[2], "Maximum coordinate in z direction")
                     .parse(in)) {
                // Only printing help returning
                return;
            }

            if (parallelism == "replicated") {
                mesh_.init_replicated();
            } else if (parallelism == "serial") {
                mesh_.init_serial();
            } else {
                mesh_.init_distributed();
            }

            auto &u_mesh = static_cast<libMesh::UnstructuredMesh &>(mesh_.raw_type());

            if (type == "file") {
            } else if (type == "sphere") {
            } else {
                libMesh::MeshTools::Generation::build_cube(
                    u_mesh,
                    b.n[0],
                    b.n[1],
                    b.n[2],
                    b.min[0],
                    b.max[0],
                    b.min[1],
                    b.max[1],
                    b.min[2],
                    b.max[2],
                    libMesh::Utility::string_to_enum<libMesh::ElemType>(elem_type),
                    gausss_lobatto_grid);
            }
            // std::shared_ptr<libMesh::MeshBase> mesh_ptr;
            // mesh_.wrap(mesh_ptr);
        }

    }  // namespace libmesh
}  // namespace utopia
