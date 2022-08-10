#include "utopia_petsc_Decompose.hpp"

#include "utopia_Instance.hpp"
#include "utopia_Logger.hpp"

#ifdef UTOPIA_WITH_METIS
#include "utopia_Metis.hpp"
#endif

#ifdef UTOPIA_WITH_PARMETIS
#include "utopia_ParMetis.hpp"
#endif

#include "utopia_petsc_CrsView.hpp"
#include "utopia_petsc_Matrix.hpp"

namespace utopia {

#ifndef UTOPIA_WITH_METIS

bool decompose(const PetscMatrix&, const int, int*)
{
    return false;
}

#else

bool decompose(const PetscMatrix& matrix, const int num_partitions, int* partitions)
{
    if (matrix.comm().size() != 1) {
        assert(false);
        Utopia::Abort("decompose only works with serial MATIJ!");
    }

    PetscCrsView mat_view = crs_view(matrix);

    idx_t nvtxs = mat_view.rows();
    idx_t ncon = 1;
    idx_t* rowptr = (idx_t*)&mat_view.row_ptr()[0];
    idx_t* colidx = (idx_t*)&mat_view.colidx()[0];
    idx_t* vwgt = nullptr;
    idx_t* vsize = nullptr;
    idx_t* adjwgt = nullptr;
    idx_t nparts = num_partitions;
    real_t* tpwgts = nullptr;
    real_t* ubvec = nullptr;
    idx_t* options = nullptr;
    idx_t objval = -1;
    idx_t* parts = partitions;

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

#ifndef UTOPIA_WITH_PARMETIS

bool parallel_decompose(const PetscMatrix&, const int, int*)
{
    return false;
}

#else

bool parallel_decompose(const PetscMatrix& matrix, const int num_partitions, int* partitions)
{
    idx_t* vtxdist = (idx_t*)&matrix.row_ranges()[0];
    idx_t ncon = 1;
    idx_t* vwgt = nullptr;
    idx_t* vsize = nullptr;
    idx_t* adjwgt = nullptr;
    idx_t nparts = num_partitions;
    real_t* tpwgts = nullptr;
    real_t* ubvec = nullptr;
    idx_t options[3] = { 0 };
    idx_t objval = -1;
    idx_t edgecut;
    idx_t* parts = partitions;
    MPI_Comm comm = matrix.comm().raw_comm();

    Mat d, o;

    const PetscInt* colmap;
    MatMPIAIJGetSeqAIJ(matrix.raw_type(), &d, &o, &colmap);

    PetscCrsView d_view(d);
    PetscCrsView o_view(o);

    std::vector<idx_t> rowptr(d_view.rows() + 1, 0);
    std::vector<idx_t> colidx(d_view.colidx().size() + o_view.colidx().size(), -1);

    auto cr = matrix.col_range();

    PetscInt local_rows = d_view.rows();
    PetscInt col_offset = cr.begin();

    auto d_rowptr = d_view.row_ptr();
    rowptr[0] = 0;
    for (PetscInt i = 0; i < local_rows; ++i) {
        rowptr[i + 1] = d_rowptr[i + 1] - d_rowptr[i];
    }

    auto o_rowptr = o_view.row_ptr();
    for (PetscInt i = 0; i < local_rows; ++i) {
        rowptr[i + 1] += o_rowptr[i + 1] - o_rowptr[i];
    }

    for (PetscInt i = 0; i < local_rows; ++i) {
        rowptr[i + 1] += rowptr[i];
    }

    for (PetscInt i = 0; i < local_rows; ++i) {
        auto d_row = d_view.row(i);
        auto o_row = o_view.row(i);

        auto begin = rowptr[i];
        for (PetscInt k = 0; k < d_row.length; ++k) {
            colidx[begin + k] = col_offset + d_row.colidx(k);
        }

        // append offsets
        begin += d_row.length;
        for (PetscInt k = 0; k < o_row.length; ++k) {
            colidx[begin + k] = colmap[o_row.colidx(k)];
        }
    }

    int ret = ParMETIS_V3_PartKway(
        vtxdist, &rowptr[0], &colidx[0], vwgt, adjwgt, 0, 0, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, parts, &comm);

    if (ret == METIS_OK) {
        return true;
    } else {
        return false;
    }
}
#endif

} // namespace utopia
