#include "utopia_sfem_SDF.hpp"

#include "matrixio_array.h"
#include "matrixio_ndarray.h"

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

            bool interpolate{true};
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
                stride[1] = nglobal[0];
                stride[2] = nglobal[0] * nglobal[1];

                in.get("interpolate", interpolate);
            }

            ~Impl() {
                if (sdf) free(sdf);
            }
        };

        SDF::SDF() : impl_(utopia::make_unique<Impl>()) {}

        SDF::~SDF() {}

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

        bool SDF::to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field) {
            auto m = (mesh_t *)mesh.raw_type();

            real_t *xnormal = (real_t *)malloc(m->nnodes * sizeof(real_t));
            real_t *ynormal = (real_t *)malloc(m->nnodes * sizeof(real_t));
            real_t *znormal = (real_t *)malloc(m->nnodes * sizeof(real_t));

            mesh.create_vector_nodal(field, 1);
            mesh.create_vector_nodal(grad_field, 3);

            auto field_view = local_view_device(field);

            if (impl_->interpolate) {
                interpolate_gap(impl_->comm.get(),
                                // Mesh
                                (ElemType)m->element_type,
                                m->nelements,
                                m->nnodes,
                                m->elements,
                                m->points,
                                // SDF
                                impl_->nlocal,
                                impl_->nglobal,
                                impl_->stride,
                                impl_->origin,
                                impl_->delta,
                                impl_->sdf,
                                // Output
                                field_view.array().begin(),
                                xnormal,
                                ynormal,
                                znormal);
            } else {
                resample_gap(impl_->comm.get(),
                             // Mesh
                             (ElemType)m->element_type,
                             m->nelements,
                             m->nnodes,
                             m->elements,
                             m->points,
                             // SDF
                             impl_->nlocal,
                             impl_->nglobal,
                             impl_->stride,
                             impl_->origin,
                             impl_->delta,
                             impl_->sdf,
                             // Output
                             field_view.array().begin(),
                             xnormal,
                             ynormal,
                             znormal);
            }

            auto grad_field_view = local_view_device(grad_field);
            for (ptrdiff_t i = 0; i < m->n_owned_nodes; i++) {
                // Convert to the vector
                grad_field_view.set(i * 3 + 0, xnormal[i]);
                grad_field_view.set(i * 3 + 1, ynormal[i]);
                grad_field_view.set(i * 3 + 2, znormal[i]);
            }

            free(xnormal);
            free(ynormal);
            free(znormal);
            return false;
        }

    }  // namespace sfem
}  // namespace utopia