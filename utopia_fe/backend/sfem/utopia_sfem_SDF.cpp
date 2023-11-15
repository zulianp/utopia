#include "utopia_sfem_SDF.hpp"

#include "matrixio_array.h"
#include "matrixio_ndarray.h"

#include "mass.h"
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
            ptrdiff_t n;

            geom_t origin[3];
            geom_t delta[3];
            geom_t *sdf{nullptr};

            Communicator comm;
            Vector weights;

            bool interpolate{false};
            bool normalize_gradient{true};

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
                sdf = (geom_t *)malloc(n * sizeof(geom_t));

                // Invalid local sizes!
                nlocal[0] = -1;
                nlocal[1] = -1;
                nlocal[2] = -1;

                if (SFEM_OK != ndarray_read(comm.get(), path.c_str(), SFEM_MPI_GEOM_T, 3, sdf, nlocal, nglobal)) {
                    utopia::err() << "Unable to read sdf file at " << path << "\n";
                    Utopia::Abort();
                }

                stride[0] = 1;
                stride[1] = nlocal[0];
                stride[2] = nlocal[0] * nlocal[1];

                // interpolate = mpi_world_size() > 1;  // FIXME!
                in.get("interpolate", interpolate);
            }

            ~Impl() {
                if (sdf) free(sdf);
            }
        };

        SDF::SDF() : impl_(utopia::make_unique<Impl>()) {}

        SDF::~SDF() {}

        bool SDF::has_weights() const { return !impl_->weights.empty(); }
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

        void SDF::clear() { impl_->weights.clear(); }

        bool SDF::interpolate_to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field) {
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

        bool SDF::project_to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field, Vector &weights) {
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

            if (mesh.comm().size() == 1) {
                resample_gap(
                    // Mesh
                    (ElemType)m->element_type,
                    m->nelements,
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
                resample_gap_local(
                    // Mesh
                    (ElemType)m->element_type,
                    m->nelements,
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

                // Utopia::Abort("...");
                assemble_lumped_mass((ElemType)m->element_type,
                                     m->nelements,
                                     m->nnodes,
                                     m->elements,
                                     m->points,
                                     weights_view.array().begin());
            }

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

        bool SDF::to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field) {
            if (interpolate()) {
                return interpolate_to_mesh(mesh, field, grad_field);
            } else {
                return project_to_mesh(mesh, field, grad_field, impl_->weights);
            }
        }

    }  // namespace sfem
}  // namespace utopia
