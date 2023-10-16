#include "utopia_sfem_SDF.hpp"

#include "matrixio_array.h"
#include "matrixio_ndarray.h"

#include "sfem_defs.h"

namespace utopia {
    namespace sfem {
        class SDF::Impl {
        public:
            ptrdiff_t nlocal[3];
            ptrdiff_t nglobal[3];
            ptrdiff_t n;

            geom_t origin[3];
            geom_t delta[3];
            geom_t *sdf{nullptr};

            Communicator comm;
        };

        SDF::SDF() {}

        SDF::~SDF() {}

        void SDF::read(Input &in) {
            // TODO
            Path path;
            in.require("path", path);

            in.require("nx", impl_->nglobal[0]);
            in.require("ny", impl_->nglobal[1]);
            in.require("nz", impl_->nglobal[2]);

            in.require("ox", impl_->origin[0]);
            in.require("oy", impl_->origin[1]);
            in.require("oz", impl_->origin[2]);

            in.require("dx", impl_->delta[0]);
            in.require("dy", impl_->delta[1]);
            in.require("dz", impl_->delta[2]);

            impl_->n = impl_->nglobal[0] * impl_->nglobal[1] * impl_->nglobal[2];
            impl_->sdf = (geom_t *)malloc(impl_->n * sizeof(geom_t));

            if (ndarray_read(
                    impl_->comm.get(), path.c_str(), SFEM_MPI_GEOM_T, 3, impl_->sdf, impl_->nlocal, impl_->nglobal)) {
                utopia::err() << "Unable to read sdf file at " << path << "\n";
                Utopia::Abort();
            }
        }

        void SDF::describe(std::ostream &os) const {
            // TODO
        }

    }  // namespace sfem
}  // namespace utopia
