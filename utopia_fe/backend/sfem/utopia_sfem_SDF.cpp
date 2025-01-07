#include "utopia_sfem_SDF.hpp"
#include <utopia_MPI.hpp>

#include "matrixio_array.h"
#include "matrixio_ndarray.h"

#include "mass.h"
#include "mesh_utils.h"
#include "node_interpolate.h"
#include "sfem_defs.h"
#include "sfem_mesh.h"
#include "sfem_resample_gap.h"

#include "utopia_sfem_Mesh.hpp"

#include "utopia_ui.hpp"

namespace utopia {
    namespace sfem {
        class SDF::Impl {
        public:
            ptrdiff_t nlocal[3];
            ptrdiff_t nglobal[3];
            ptrdiff_t stride[3];
            ptrdiff_t n{0};

            geom_t origin[3];
            geom_t delta[3];
            geom_t *sdf{nullptr};

            Communicator comm;
            Vector weights;

            bool interpolate{false};
            bool normalize_gradient{true};
            bool has_weights{false};

            // Sharp feature stuff
            bool geometry_aware_resampling{false};
            bool super_imposed{false};
            geom_t angle_threshold{0.15};

            void read(Input &in) {
                Path path;
                in.require("path", path);

                in.require("nx", nglobal[0]);
                in.require("ny", nglobal[1]);
                in.require("nz", nglobal[2]);

                in.require("ox", origin[0]);
                in.require("oy", origin[1]);
                in.require("oz", origin[2]);

                in.require("dx", delta[0]);
                in.require("dy", delta[1]);
                in.require("dz", delta[2]);

                n = nglobal[0] * nglobal[1] * nglobal[2];
                sdf = nullptr;

                // Invalid local sizes!
                nlocal[0] = -1;
                nlocal[1] = -1;
                nlocal[2] = -1;

                if (SFEM_OK != ndarray_create_from_file(
                                   comm.get(), path.c_str(), SFEM_MPI_GEOM_T, 3, (void **)&sdf, nlocal, nglobal)) {
                    utopia::err() << "Unable to read sdf file at " << path << "\n";
                    Utopia::Abort("IO error!");
                }

                stride[0] = 1;
                stride[1] = nlocal[0];
                stride[2] = nlocal[0] * nlocal[1];

                // interpolate = mpi_world_size() > 1;  // FIXME!
                in.get("interpolate", interpolate);
                in.get("angle_threshold", angle_threshold);
                in.get("geometry_aware_resampling", geometry_aware_resampling);
                in.get("super_imposed", super_imposed);

                bool verbose = false;
                in.get("verbose", verbose);

                if (verbose && !mpi_world_rank()) {
                    utopia::out() << "utopia::sfem::SDF: interpolate = " << interpolate << '\n';
                }
            }

            ~Impl() {
                if (sdf) free(sdf);
            }
        };

        SDF::SDF() : impl_(utopia::make_unique<Impl>()) {}

        SDF::~SDF() {}

        bool SDF::has_weights() const { return impl_->has_weights; }
        const SDF::Vector &SDF::weights() const { return impl_->weights; }

        void SDF::read_from_file(const Path &path) {
            auto input = open_istream(path);
            impl_->read(*input);
        }

        void SDF::read(Input &in) {
            Path path;
            in.require("path", path);

            if (path.extension() == "yaml" || path.extension() == "yml" || path.extension() == "json") {
                read_from_file(path);
                return;
            } else {
                impl_->read(in);
            }
        }

        void SDF::describe(std::ostream &os) const {
            // TODO
        }

        bool SDF::interpolate() const { return impl_->interpolate; }

        void SDF::clear() {
            impl_->weights.clear();
            impl_->has_weights = false;
        }

        bool SDF::interpolate_to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field) {
            UTOPIA_TRACE_SCOPE("SDF::interpolate_to_mesh");
            clear();

            auto m = (mesh_t *)mesh.raw_type();

            real_t *xnormal = (real_t *)malloc(m->nnodes * sizeof(real_t));
            real_t *ynormal = (real_t *)malloc(m->nnodes * sizeof(real_t));
            real_t *znormal = (real_t *)malloc(m->nnodes * sizeof(real_t));

            mesh.create_vector_nodal(field, 1);
            mesh.create_vector_nodal(grad_field, 3);

            auto field_view = local_view_device(field);

            geom_t *actual_sdf = nullptr;
            geom_t *psdf = nullptr;

            ptrdiff_t nlocal[3] = {impl_->nlocal[0], impl_->nlocal[1], impl_->nlocal[2]};
            geom_t origin[3] = {impl_->origin[0], impl_->origin[1], impl_->origin[2]};

            if (mesh.comm().size() > 1) {
                sdf_view(mesh.comm().get(),
                         m->nnodes,
                         m->points[2],
                         impl_->nlocal,
                         impl_->nglobal,
                         impl_->stride,
                         impl_->origin,
                         impl_->delta,
                         impl_->sdf,
                         &psdf,
                         &nlocal[2],
                         &origin[2]);

                actual_sdf = psdf;
            } else {
                actual_sdf = impl_->sdf;
            }

            interpolate_gap(m->n_owned_nodes,
                            m->points,
                            // SDF
                            nlocal,
                            impl_->stride,
                            origin,
                            impl_->delta,
                            actual_sdf,
                            // Output
                            field_view.array().begin(),
                            xnormal,
                            ynormal,
                            znormal);

            auto grad_field_view = local_view_device(grad_field);
            for (ptrdiff_t i = 0; i < m->n_owned_nodes; i++) {
                // Convert to the vector
                grad_field_view.set(i * 3 + 0, xnormal[i]);
                grad_field_view.set(i * 3 + 1, ynormal[i]);
                grad_field_view.set(i * 3 + 2, znormal[i]);
            }

            if (psdf) {
                free(psdf);
            }

            free(xnormal);
            free(ynormal);
            free(znormal);
            return false;
        }

        bool SDF::project_to_mesh_with_sharp_features(const Mesh &mesh,
                                                      Vector &field,
                                                      Vector &grad_field,
                                                      Vector &weights) {
            UTOPIA_TRACE_SCOPE("SDF::project_to_mesh_with_sharp_features");

            clear();
            impl_->has_weights = true;

            auto m = (mesh_t *)mesh.raw_type();

            real_t *xnormal = (real_t *)calloc(m->nnodes, sizeof(real_t));
            real_t *ynormal = (real_t *)calloc(m->nnodes, sizeof(real_t));
            real_t *znormal = (real_t *)calloc(m->nnodes, sizeof(real_t));

            mesh.create_vector_nodal(field, 1);
            mesh.create_vector_nodal(grad_field, 3);

            field.set(0);
            grad_field.set(0);

            auto field_view = local_view_device(field);

            geom_t *actual_sdf = nullptr;
            geom_t *psdf = nullptr;

            ptrdiff_t nlocal[3] = {impl_->nlocal[0], impl_->nlocal[1], impl_->nlocal[2]};
            geom_t origin[3] = {impl_->origin[0], impl_->origin[1], impl_->origin[2]};

            if (mesh.comm().size() > 1) {
                sdf_view(mesh.comm().get(),
                         m->nnodes,
                         m->points[2],
                         impl_->nlocal,
                         impl_->nglobal,
                         impl_->stride,
                         impl_->origin,
                         impl_->delta,
                         impl_->sdf,
                         &psdf,
                         &nlocal[2],
                         &origin[2]);

                actual_sdf = psdf;
            } else {
                actual_sdf = impl_->sdf;
            }

            auto edges_ptr = mesh.sharp_edges(impl_->angle_threshold);
            auto edge_elements = edges_ptr->elements().reinterpret<idx_t *>();

            weights.zeros(layout(field));
            auto weights_view = local_view_device(weights);

            {
                // Edges
                resample_gap_local(
                    // Mesh
                    shell_type((ElemType)edges_ptr->element_type()),
                    edges_ptr->n_elements(),
                    m->nnodes,
                    edge_elements.begin(),
                    m->points,
                    // SDF
                    nlocal,
                    impl_->stride,
                    origin,
                    impl_->delta,
                    actual_sdf,
                    // Output
                    field_view.array().begin(),
                    xnormal,
                    ynormal,
                    znormal);

                assemble_lumped_mass(shell_type((ElemType)edges_ptr->element_type()),
                                     edges_ptr->n_elements(),
                                     m->nnodes,
                                     edge_elements.begin(),
                                     m->points,
                                     weights_view.array().begin());
            }

            {
                // Faces
                if (impl_->super_imposed) {
                    // Faces (all)
                    resample_gap_local(
                        // Mesh
                        shell_type((ElemType)m->element_type),
                        m->n_owned_elements,
                        m->nnodes,
                        m->elements,
                        m->points,
                        // SDF
                        nlocal,
                        impl_->stride,
                        origin,
                        impl_->delta,
                        actual_sdf,
                        // Output
                        field_view.array().begin(),
                        xnormal,
                        ynormal,
                        znormal);

                    assemble_lumped_mass(shell_type((ElemType)m->element_type),
                                         m->n_owned_elements,
                                         m->nnodes,
                                         m->elements,
                                         m->points,
                                         weights_view.array().begin());
                } else {
                    // Faces (disconnected)
                    auto faces_ptr = mesh.disconnected_faces_from_sharp_edges(*edges_ptr);

                    if (faces_ptr) {
                        auto face_elements = faces_ptr->elements().reinterpret<idx_t *>();

                        resample_gap_local(
                            // Mesh
                            shell_type((ElemType)faces_ptr->element_type()),
                            faces_ptr->n_elements(),
                            m->nnodes,
                            face_elements.begin(),
                            m->points,
                            // SDF
                            nlocal,
                            impl_->stride,
                            origin,
                            impl_->delta,
                            actual_sdf,
                            // Output
                            field_view.array().begin(),
                            xnormal,
                            ynormal,
                            znormal);

                        assemble_lumped_mass(shell_type((ElemType)faces_ptr->element_type()),
                                             faces_ptr->n_elements(),
                                             m->nnodes,
                                             face_elements.begin(),
                                             m->points,
                                             weights_view.array().begin());
                    }
                }
            }

            {
                // Corners
                auto corners_ptr = mesh.separate_corners_from_sharp_edges(*edges_ptr, !impl_->super_imposed);

                if (corners_ptr) {
                    auto corner_elements = corners_ptr->elements().reinterpret<idx_t *>();

                    real_t *p_g = (real_t *)calloc(corners_ptr->n_elements(), sizeof(real_t));
                    real_t *p_xnormal = (real_t *)calloc(corners_ptr->n_elements(), sizeof(real_t));
                    real_t *p_ynormal = (real_t *)calloc(corners_ptr->n_elements(), sizeof(real_t));
                    real_t *p_znormal = (real_t *)calloc(corners_ptr->n_elements(), sizeof(real_t));

                    geom_t **corner_points = allocate_points(m->spatial_dim, corners_ptr->n_elements());
                    select_points(
                        m->spatial_dim, corners_ptr->n_elements(), corner_elements[0], m->points, corner_points);

                    interpolate_gap(
                        // Mesh
                        corners_ptr->n_elements(),
                        corner_points,
                        // SDF
                        nlocal,
                        impl_->stride,
                        origin,
                        impl_->delta,
                        actual_sdf,
                        // Output
                        p_g,
                        p_xnormal,
                        p_ynormal,
                        p_znormal);

                    // Add to complete array

                    auto mass_vector = weights_view.array();
                    auto g = field_view.array();

                    for (ptrdiff_t i = 0; i < corners_ptr->n_elements(); i++) {
                        g[corner_elements[0][i]] += p_g[i];
                        xnormal[corner_elements[0][i]] += p_xnormal[i];
                        ynormal[corner_elements[0][i]] += p_ynormal[i];
                        znormal[corner_elements[0][i]] += p_znormal[i];
                        mass_vector[corner_elements[0][i]] += 1;
                    }

                    free_points(m->spatial_dim, corner_points);
                    free(p_g);
                    free(p_xnormal);
                    free(p_ynormal);
                    free(p_znormal);
                }
            }

            auto grad_field_view = local_view_device(grad_field);
            // for (ptrdiff_t i = 0; i < m->n_owned_nodes; i++) {
            for (ptrdiff_t i = 0; i < m->nnodes; i++) {
                // Convert to the vector
                grad_field_view.set(i * 3 + 0, xnormal[i]);
                grad_field_view.set(i * 3 + 1, ynormal[i]);
                grad_field_view.set(i * 3 + 2, znormal[i]);
            }

            if (psdf) {
                free(psdf);
            }

            free(xnormal);
            free(ynormal);
            free(znormal);
            return false;
        }

        bool SDF::project_to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field, Vector &weights) {
            UTOPIA_TRACE_SCOPE("SDF::project_to_mesh");

            if (impl_->geometry_aware_resampling) {
                return project_to_mesh_with_sharp_features(mesh, field, grad_field, weights);
            }

            clear();

            auto m = (mesh_t *)mesh.raw_type();

            real_t *xnormal = (real_t *)calloc(m->nnodes, sizeof(real_t));
            real_t *ynormal = (real_t *)calloc(m->nnodes, sizeof(real_t));
            real_t *znormal = (real_t *)calloc(m->nnodes, sizeof(real_t));

            mesh.create_vector_nodal(field, 1);
            mesh.create_vector_nodal(grad_field, 3);

            field.set(0);
            grad_field.set(0);

            auto field_view = local_view_device(field);

            geom_t *actual_sdf = nullptr;
            geom_t *psdf = nullptr;

            ptrdiff_t nlocal[3] = {impl_->nlocal[0], impl_->nlocal[1], impl_->nlocal[2]};
            geom_t origin[3] = {impl_->origin[0], impl_->origin[1], impl_->origin[2]};

            if (mesh.comm().size() > 1) {
                sdf_view(mesh.comm().get(),
                         m->nnodes,
                         m->points[2],
                         impl_->nlocal,
                         impl_->nglobal,
                         impl_->stride,
                         impl_->origin,
                         impl_->delta,
                         impl_->sdf,
                         &psdf,
                         &nlocal[2],
                         &origin[2]);

                actual_sdf = psdf;
            } else {
                actual_sdf = impl_->sdf;
            }

            if (mesh.comm().size() == 1) {
                resample_gap(
                    // Mesh
                    (ElemType)m->element_type,
                    m->n_owned_elements,
                    m->nnodes,
                    m->elements,
                    m->points,
                    // SDF
                    nlocal,
                    impl_->stride,
                    origin,
                    impl_->delta,
                    actual_sdf,
                    // Output
                    field_view.array().begin(),
                    xnormal,
                    ynormal,
                    znormal);

            } else {
                impl_->has_weights = true;

                resample_gap_local(
                    // Mesh
                    (ElemType)m->element_type,
                    m->n_owned_elements,
                    m->nnodes,
                    m->elements,
                    m->points,
                    // SDF
                    nlocal,
                    impl_->stride,
                    origin,
                    impl_->delta,
                    actual_sdf,
                    // Output
                    field_view.array().begin(),
                    xnormal,
                    ynormal,
                    znormal);

                weights.zeros(layout(field));

                auto weights_view = local_view_device(weights);

                assemble_lumped_mass(shell_type((ElemType)m->element_type),
                                     m->n_owned_elements,
                                     m->nnodes,
                                     m->elements,
                                     m->points,
                                     weights_view.array().begin());
            }

            if (m->n_owned_elements) {
                auto grad_field_view = local_view_device(grad_field);
                // for (ptrdiff_t i = 0; i < m->n_owned_nodes; i++) {
                for (ptrdiff_t i = 0; i < m->nnodes; i++) {
                    // Convert to the vector
                    grad_field_view.set(i * 3 + 0, xnormal[i]);
                    grad_field_view.set(i * 3 + 1, ynormal[i]);
                    grad_field_view.set(i * 3 + 2, znormal[i]);
                }
            }

            if (psdf) {
                free(psdf);
            }

            free(xnormal);
            free(ynormal);
            free(znormal);
            return false;
        }

        bool SDF::to_mesh(const Mesh &mesh, Vector &field) {
            Vector dump;
            return to_mesh(mesh, field, dump);
        }

        bool SDF::to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field) {
            if (interpolate()) {
                return interpolate_to_mesh(mesh, field, grad_field);
            } else {
                return project_to_mesh(mesh, field, grad_field, impl_->weights);
            }
        }

    }  // namespace sfem
}  // namespace utopia
