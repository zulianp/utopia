#include "utopia_libmesh_MeshInitializer.hpp"
#include "utopia_Options.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_libmesh_Mesh.hpp"

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

            void make_square() {
                n[2] = 0;
                min[2] = max[2] = 0.0;
            }

            void make_line() {
                n[1] = n[2] = 0;
                min[1] = max[1] = 0.0;
                min[2] = max[2] = 0.0;
            }

            void zero() { default_cube(0, 0.0, 0.0); }

            int n[max_dim];
            libMesh::Real min[max_dim], max[max_dim];
        };

        class MeshInitializer::Scale : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options().add_option("scale", scale, "Uniform scaling the mesh.").parse(in)) {
                    return;
                }
            }

            void describe(std::ostream &os) const override { os << "Scale"; }

            void apply(libMesh::UnstructuredMesh &mesh) {
                if (1. == scale) return;

                for (auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {
                    for (int i = 0; i < LIBMESH_DIM; ++i) {
                        (**it)(i) *= scale;
                    }
                }
            }

            libMesh::Real scale{1.0};
        };

        class MeshInitializer::Shift : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options()
                         .add_option("shift_x", shift[0], "Shift in x direction.")
                         .add_option("shift_y", shift[1], "Shift in y direction.")
                         .add_option("shift_z", shift[2], "Shift in z direction.")
                         .parse(in)) {
                    return;
                }
            }

            void describe(std::ostream &os) const override { os << "Shift"; }

            void apply(libMesh::UnstructuredMesh &mesh) {
                if (0. == shift[0] && 0. == shift[1] && 0. == shift[2]) return;

                for (auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {
                    for (int i = 0; i < LIBMESH_DIM; ++i) {
                        (**it)(i) += shift[i];
                    }
                }
            }

            libMesh::Real shift[3] = {0.0, 0.0, 0.0};
        };

        MeshInitializer::MeshInitializer(Mesh &mesh) : mesh_(mesh) {}

        void MeshInitializer::read(Input &in) {
            std::string parallelism = "distributed";
            std::string type = "cube";
            std::string elem_type = "HEX8";
            std::string path = "";
            int n = 10;
            libMesh::Real min = 0.0, max = 1.0;
            bool gauss_lobatto_grid = false;
            libMesh::Real radius = 1.0;
            int refinements = 0;
            bool all_tri = false;
            bool all_second_order = false;
            bool full_order = true;

            Build b;
            b.default_cube(n, min, max);

            Shift shift;
            Scale scale;

            if (!Options()
                     .add_option(
                         "parallelism", parallelism, "Parallelism avaialble in libMesh: distributed|replicated|serial.")
                     .add_option("type", type, "Type of the mesh: line|square|cube|sphere|file.")
                     .add_option("path", path, "Path of the mesh when type=file.")
                     .add_option("elem_type", elem_type, "Mesh element type. See libMesh ElemType enum.")
                     // .add_option("n", n, "Number of nodes for all dimensions.")
                     // .add_option("min", min, "Minumum coordinate for all directions.")
                     // .add_option("max", max, "Maximum coordinate for all directions.")
                     .add_option("nx", b.n[0], "Number of nodes in x direction.")
                     .add_option("ny", b.n[1], "Number of nodes in y direction.")
                     .add_option("nz", b.n[2], "Number of nodes in z direction.")
                     .add_option("min_x", b.min[0], "Minumum coordinate in x direction.")
                     .add_option("min_y", b.min[1], "Minumum coordinate in y direction.")
                     .add_option("min_z", b.min[2], "Minumum coordinate in z direction.")
                     .add_option("max_x", b.max[0], "Maximum coordinate in x direction.")
                     .add_option("max_y", b.max[1], "Maximum coordinate in y direction.")
                     .add_option("max_z", b.max[2], "Maximum coordinate in z direction.")
                     .add_option("gauss_lobatto_grid", gauss_lobatto_grid, "Use Gauss-Lobatto grid.")
                     .add_option("radius", radius, "Sphere radius.")
                     .add_option("refinements", refinements, "Number of mesh refinements.")
                     .add_option("all_tri", all_tri, "Convert to triangle mesh.")
                     .add_option("all_second_order", all_second_order, "Convert to second order mesh when type=file.")
                     .add_option("full_order",
                                 full_order,
                                 "Convert to full second order mesh when all_second_order=true and type=file.")
                     .add_option(
                         "shift_x", shift.shift[0], "Applies a shift transform to the mesh in the x coordinate.")
                     .add_option(
                         "shift_y", shift.shift[1], "Applies a shift transform to the mesh in the y coordinate.")
                     .add_option(
                         "shift_z", shift.shift[2], "Applies a shift transform to the mesh in the z coordinate.")
                     .add_option("scale", scale.scale, "Applies a scale transform to the mesh.")
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
                if (path.empty()) {
                    err() << "[Error] path for file is undefined\n";
                    Utopia::Abort();
                } else {
                    u_mesh.read(path);
                }
            } else if (type == "sphere") {
                libMesh::MeshTools::Generation::build_sphere(
                    u_mesh,
                    radius,
                    std::max(1, refinements),
                    libMesh::Utility::string_to_enum<libMesh::ElemType>(elem_type)
                    // const unsigned int    n_smooth = 2,
                    // const bool   flat = true
                );
            } else {
                if (type == "line") {
                    b.make_line();
                } else if (type == "square") {
                    b.make_square();
                }

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
                    gauss_lobatto_grid);
            }

            if (all_tri) {
                libMesh::MeshTools::Modification::all_tri(u_mesh);
            }

            if (type == "file" && all_second_order) {
                u_mesh.all_second_order(full_order);
            }

            scale.apply(u_mesh);
            shift.apply(u_mesh);
        }

    }  // namespace libmesh
}  // namespace utopia
