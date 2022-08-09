#include "utopia_petsc_Decompose.hpp"

#include "utopia_Instance.hpp"
#include "utopia_Logger.hpp"

#include "utopia_Metis.hpp"

#include "utopia_petsc_CrsView.hpp"
#include "utopia_petsc_Matrix.hpp"

namespace utopia {

#ifndef UTOPIA_WITH_METIS

    bool decompose(const PetscMatrix &, const int, int *) { return false; }

#else

    bool decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions) {
        if (matrix.comm().size() != 1) {
            assert(false);
            Utopia::Abort("decompose only works with serial MATIJ!");
        }

        PetscCrsView mat_view = crs_view(matrix);

        idx_t nvtxs = mat_view.rows();
        idx_t ncon = 1;
        idx_t *rowptr = (idx_t *)&mat_view.row_ptr()[0];
        idx_t *colidx = (idx_t *)&mat_view.colidx()[0];
        idx_t *vwgt = nullptr;
        idx_t *vsize = nullptr;
        idx_t *adjwgt = nullptr;
        idx_t nparts = num_partitions;
        real_t *tpwgts = nullptr;
        real_t *ubvec = nullptr;
        idx_t *options = nullptr;
        idx_t objval = -1;
        idx_t *parts = partitions;

        int ret = METIS_PartGraphKway(
            &nvtxs, &ncon, rowptr, colidx, vwgt, vsize, adjwgt, &nparts, tpwgts, ubvec, options, &objval, parts);

        if (ret == METIS_OK) {
            return true;
        } else {
            std::string message;
            switch (ret) {
                case METIS_ERROR_INPUT: {
                    message = "METIS_ERROR_INPUT";
                    break;
                }
                case METIS_ERROR_MEMORY: {
                    message = "METIS_ERROR_MEMORY";
                    break;
                }
                case METIS_ERROR: {
                    message = "METIS_ERROR";
                    break;
                }
                default: {
                    message = "Unknow error";
                }
            }

            utopia::err() << "Metis return error code " << ret << " (" << message << ")\n";
            return false;
        }
    }

#endif
}  // namespace utopia
