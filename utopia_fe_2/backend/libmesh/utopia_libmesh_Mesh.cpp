#include "utopia_libmesh_Mesh.hpp"

// utopia
#include "utopia_make_unique.hpp"

// utopia fe

// libmesh
#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/string_to_enum.h"

namespace utopia {

    Mesh<libMesh::UnstructuredMesh>::Mesh(const Communicator &comm)
        : comm_(std::make_shared<libMesh::Parallel::Communicator>(comm.raw_comm())),
          impl_(utopia::make_unique<libMesh::DistributedMesh>(*comm_)) {}

    Mesh<libMesh::UnstructuredMesh>::~Mesh() {}

    void Mesh<libMesh::UnstructuredMesh>::describe(std::ostream &os) const {
        os << impl_->mesh_dimension() << std::endl;
    }

    void Mesh<libMesh::UnstructuredMesh>::read(Input &in) {
        std::string mesh_type = "square";
        std::string path = "";

        int refinements = 0;

        Scalar span[3] = {0., 0., 0.};

        Scalar min_coords[3] = {0., 0., 0.};
        Scalar max_coords[3] = {1., 1., 1.};

        int n[3] = {5, 5, 5};
        int order = 1;

        Scalar scale = 1.;
        Scalar shift[3] = {0., 0., 0.};

        bool full_order = true;

        std::string elem_type = "quad";

        in.get("type", mesh_type);
        in.get("elem_type", elem_type);
        in.get("order", order);
        in.get("full_order", full_order);
        in.get("path", path);

        in.get("refinements", refinements);

        in.get("span_x", span[0]);
        in.get("span_y", span[1]);
        in.get("span_z", span[2]);

        in.get("min_x", min_coords[0]);
        in.get("min_y", min_coords[1]);
        in.get("min_z", min_coords[2]);

        in.get("max_x", max_coords[0]);
        in.get("max_y", max_coords[1]);
        in.get("max_z", max_coords[2]);

        in.get("n_x", n[0]);
        in.get("n_y", n[1]);
        in.get("n_z", n[2]);

        in.get("scale", scale);
        in.get("shift_x", shift[0]);
        in.get("shift_y", shift[1]);
        in.get("shift_z", shift[2]);

        libMesh::ElemType lm_elem_type = libMesh::INVALID_ELEM;

        if (!elem_type.empty()) {
            lm_elem_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(elem_type);
        }

        if (mesh_type == "file") {
            impl_->read(path);
        } else if (mesh_type == "line") {
            libMesh::MeshTools::Generation::build_line(*impl_, n[0], min_coords[0], max_coords[0], lm_elem_type);
        } else if (mesh_type == "square") {
            libMesh::MeshTools::Generation::build_square(
                *impl_, n[0], n[1], min_coords[0], max_coords[0], min_coords[1], max_coords[1], lm_elem_type);
        } else if (mesh_type == "cube") {
            libMesh::MeshTools::Generation::build_cube(*impl_,
                                                       n[0],
                                                       n[1],
                                                       n[2],
                                                       min_coords[0],
                                                       max_coords[0],
                                                       min_coords[1],
                                                       max_coords[1],
                                                       min_coords[2],
                                                       max_coords[2],
                                                       lm_elem_type);
        } else if (mesh_type == "sphere") {
            Scalar radius = 1.;
            int sphere_refine = 2;

            in.get("radius", radius);
            in.get("sphere_refine", sphere_refine);

            libMesh::MeshTools::Generation::build_sphere(*impl_,
                                                         radius,
                                                         sphere_refine,  // const unsigned int nr = 2,
                                                         lm_elem_type
                                                         // const unsigned int    n_smooth = 2,
                                                         // const bool   flat = true
            );
        }
        // else if (mesh_type == "aabb") {
        //     libMesh::DistributedMesh temp_mesh(impl_->comm());
        //     temp_mesh.read(path);

        //     auto bb = bounding_box(temp_mesh);

        //     if (temp_mesh.spatial_dimension() == 3) {
        //         libMesh::MeshTools::Generation::build_cube(*impl_,
        //                                                    n[0],
        //                                                    n[1],
        //                                                    n[2],
        //                                                    bb.min()(0) - span[0],
        //                                                    bb.max()(0) + span[0],
        //                                                    bb.min()(1) - span[1],
        //                                                    bb.max()(1) + span[1],
        //                                                    bb.min()(2) - span[2],
        //                                                    bb.max()(2) + span[2],
        //                                                    get_type(elem_type, order, 3, full_order));

        //     } else {
        //         libMesh::MeshTools::Generation::build_square(*impl_,
        //                                                      n[0],
        //                                                      n[1],
        //                                                      bb.min()(0) - span[0],
        //                                                      bb.max()(0) + span[0],
        //                                                      bb.min()(1) - span[1],
        //                                                      bb.max()(1) + span[1],
        //                                                      get_type(elem_type, order, 2, full_order));
        //     }
        // }

        // int block_override = -1;
        // in.get("block-override", block_override);

        // override_block(block_override, *impl_);

        // // build_extrusion (UnstructuredMesh &mesh, const UnstructuredMesh &cross_section, const unsigned int nz,
        // // RealVectorValue extrusion_vector, QueryElemSubdomainIDBase *elem_subdomain=libmesh_nullptr)

        // scale_mesh(scale, *impl_);
        // shift_mesh(shift, *impl_);

        // SideSetAssignment<libMesh::UnstructuredMesh> ssa;
        // in.get("side-set-assignement", ssa);
        // ssa.apply(*impl_);

        // refine(refinements, *impl_);

        // bool convert_to_triangles = false;
        // in.get("convert-to-triangles", convert_to_triangles);

        // if (convert_to_triangles) {
        //     libMesh::MeshTools::Modification::all_tri(*impl_);
        // }

        // if (mesh_type == "file" && order == 2) {
        //     impl_->all_second_order(full_order);
        // }

        // {
        //     // single morph
        //     UIMorph<libMesh::DistributedMesh> morph;
        //     in.get("morph", morph);

        //     if (morph.is_valid()) {
        //         morph.apply(*impl_);
        //     }
        // }

        // // in.get("improve-smoothness", [this](Input &is) {
        // //     MeshParamSmoother smoother;
        // //     smoother.read(is);
        // //     smoother.apply(*impl_);
        // // });

        // // multi morph
        // in.get("morphs", [this](Input &is) {
        //     in.get_all([this](Input &is) {
        //         UIMorph<libMesh::DistributedMesh> morph;
        //         morph.read(is);

        //         if (morph.is_valid()) {
        //             morph.apply(*impl_);
        //         }
        //     });
        // });

        // bool must_refine_0 = false;
        // in.get("refine-0", must_refine_0);

        // bool must_random_refine = false;
        // in.get("random-refine", must_random_refine);

        // int n_random_refinements = 1;
        // in.get("random-refinements", n_random_refinements);

        // if (must_refine_0) {
        //     refine_0(*impl_);
        // }

        // if (must_random_refine) {
        //     random_refine(*impl_, n_random_refinements);
        // }

        // bool must_refine_at_intersection = false;
        // in.get("refine-at-intersection", must_refine_at_intersection);

        // if (must_refine_at_intersection) {
        //     std::string intersecting_mesh_path;
        //     in.get("intersecting-mesh-path", intersecting_mesh_path);

        //     int n_refinements_at_intersection = 1;
        //     in.get("n-refinements-at-intersection", n_refinements_at_intersection);

        //     auto i_mesh = std::make_shared<libMesh::DistributedMesh>(impl_->comm());

        //     i_mesh->read(intersecting_mesh_path);

        //     refine_at_intersection(i_mesh, libMesh::Order(1), mesh_, n_refinements_at_intersection, false);
        // }
    }

}  // namespace utopia

// void refine_at_intersection(const std::shared_ptr<libMesh::UnstructuredMesh> &fracture_network,
//                             const libMesh::Order &elem_order,
//                             const std::shared_ptr<libMesh::UnstructuredMesh> &mesh,
//                             const int refinement_loops = 1,
//                             const bool use_interpolation = false);

// inline void refine_0(libMesh::UnstructuredMesh &mesh) {
//     libMesh::MeshRefinement mesh_refinement(mesh);
//     auto e_it = elements_begin(mesh);

//     if (e_it != elements_end(mesh)) {
//         (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
//     }

//     mesh_refinement.make_flags_parallel_consistent();
//     mesh_refinement.refine_elements();
//     mesh_refinement.test_level_one(true);

//     mesh.prepare_for_use();
// }

// inline void random_refine(libMesh::UnstructuredMesh &mesh, const int refinement_loops) {
//     libMesh::MeshRefinement mesh_refinement(mesh);

//     for (int i = 0; i < refinement_loops; ++i) {
//         mesh_refinement.clean_refinement_flags();

//         int idx = 0;
//         for (auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it, ++idx) {
//             auto val = idx % 2 == 1;
//             if (val) {
//                 (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
//             }
//         }

//         mesh_refinement.make_flags_parallel_consistent();
//         mesh_refinement.refine_elements();
//         mesh_refinement.test_level_one(true);
//     }

//     // mesh_refinement.clean_refinement_flags()
// }

// static void refine(const int n_refs, libMesh::UnstructuredMesh &mesh) {
//     if (n_refs <= 0) return;

//     libMesh::MeshRefinement mesh_refinement(mesh);
//     mesh_refinement.make_flags_parallel_consistent();
//     mesh_refinement.uniformly_refine(n_refs);
// }

// static void scale_mesh(const Scalar &scale_factor, libMesh::UnstructuredMesh &mesh) {
//     if (scale_factor == 1.) return;
//     assert(scale_factor > 0.);

//     for (auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {
//         for (int i = 0; i < LIBMESH_DIM; ++i) {
//             (**it)(i) *= scale_factor;
//         }
//     }
// }

// static void override_block(const int block, libMesh::UnstructuredMesh &mesh) {
//     if (block <= 0) return;

//     for (auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
//         (*it)->subdomain_id() = block;
//     }
// }

// static void shift_mesh(const Scalar t[3], libMesh::UnstructuredMesh &mesh) {
//     if (0. == t[0] && 0. == t[1] && 0. == t[2]) return;

//     for (auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {
//         for (int i = 0; i < LIBMESH_DIM; ++i) {
//             (**it)(i) += t[i];
//         }
//     }
// }
