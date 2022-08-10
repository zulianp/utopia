#include "utopia_petsc_Decompose.hpp"

#include "utopia_Communicator.hpp"
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

    bool decompose(const PetscMatrix &, const int, int *) { return false; }

#else

    bool decompose(const PetscMatrix& matrix, const int num_partitions, int* partitions) {
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

    bool parallel_decompose(const PetscMatrix &, const int, int *) { return false; }

#else

    bool parallel_decompose(const PetscMatrix& matrix, const int num_partitions, int* partitions) {
        idx_t* vtxdist = (idx_t*)&matrix.row_ranges()[0];
        idx_t ncon = 1;
        idx_t* vwgt = nullptr;
        idx_t* vsize = nullptr;
        idx_t* adjwgt = nullptr;
        idx_t wgtflag = 0;
        idx_t numflag = 0;
        idx_t nparts = num_partitions;

        real_t ubvec[1] = {1.05};
        idx_t options[3] = {0};
        idx_t objval = -1;
        idx_t edgecut;
        idx_t* parts = partitions;
        MPI_Comm comm = matrix.comm().raw_comm();
        int comm_size = matrix.comm().size();

        Mat d, o;

        const PetscInt* colmap;
        MatMPIAIJGetSeqAIJ(matrix.raw_type(), &d, &o, &colmap);

        PetscCrsView d_view(d);
        PetscCrsView o_view(o);

        std::vector<idx_t> rowptr(d_view.rows() + 1, 0);
        std::vector<idx_t> colidx(d_view.colidx().size() + o_view.colidx().size(), -1);
        std::vector<real_t> tpwgts(num_partitions, 1. / num_partitions);
        std::vector<idx_t> actual_vwgts(d_view.rows(), 0);

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

        // Copy weights
        for (PetscInt i = 0; i < local_rows; ++i) {
            actual_vwgts[i] = rowptr[i + 1];
        }
        vwgt = &actual_vwgts[0];
        wgtflag = 2;

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

        if (false) {
            std::stringstream ss;

            for (auto r : rowptr) {
                ss << r << " ";
            }

            ss << "\n";

            for (auto c : colidx) {
                ss << c << " ";
            }

            ss << "\n";

            real_t sumw = 0;

            for (auto w : tpwgts) {
                sumw += w;
            }

            ss << "sumw = " << sumw << "\n";

            matrix.comm().synched_print(ss.str());
        }

        int ret = ParMETIS_V3_PartKway(vtxdist,     // 0
                                       &rowptr[0],  // 1
                                       &colidx[0],  // 2
                                       vwgt,        // 3
                                       adjwgt,      // 4
                                       &wgtflag,    // 5
                                       &numflag,    // 6
                                       &ncon,       // 7
                                       &nparts,     // 8
                                       &tpwgts[0],  // 9
                                       ubvec,       // 10
                                       options,     // 11
                                       &edgecut,    // 12
                                       parts,       // 13
                                       &comm);      // 14

        if (ret == METIS_OK) {
            return true;
        } else {
            return false;
        }
    }
#endif

    // Remote data is arranged locally following rank ordering
    bool partitions_to_permutations(const PetscMatrix &matrix,
                                    const int *partitions,
                                    Traits<PetscMatrix>::SizeType *index) {
        auto &&comm = matrix.comm();
        const int comm_size = comm.size();

        auto rrs = matrix.row_ranges();

        std::vector<PetscInt> send_row_count(comm_size, 0);
        std::vector<PetscInt> recv_row_count(comm_size, 0);

        const int local_rows = matrix.local_rows();

        int incoming = 0;

        for (int i = 0; i < local_rows; ++i) {
            send_row_count[partitions[i]]++;
        }

        // exchange send/receive information
        MPI_Alltoall(send_row_count.data(),
                     1,
                     utopia::MPIType<PetscInt>::value(),
                     recv_row_count.data(),
                     1,
                     utopia::MPIType<PetscInt>::value(),
                     comm.raw_comm());

        for (int r = 0; r < comm_size; ++r) {
            // Rows that are kept local are also counted!
            incoming += recv_row_count[r];
        }

        int proc_offset = 0;
        comm.exscan(&incoming, &proc_offset, 1, MPI_SUM);

        std::vector<PetscInt> index_buff(comm_size, 0);
        std::vector<PetscInt> index_offset(comm_size, 0);

        index_buff[0] = proc_offset;

        for (int r = 0; r < comm_size - 1; ++r) {
            index_buff[r + 1] = recv_row_count[r] + index_buff[r];
        }

        MPI_Alltoall(index_buff.data(),
                     1,
                     utopia::MPIType<PetscInt>::value(),
                     index_offset.data(),
                     1,
                     utopia::MPIType<PetscInt>::value(),
                     comm.raw_comm());

        for (int i = 0; i < local_rows; ++i) {
            index[i] = index_offset[partitions[i]]++;
        }

        return true;
    }

    bool rebalance(const PetscMatrix &in,
                   PetscMatrix &out,
                   Traits<PetscMatrix>::IndexArray &partitioning,
                   Traits<PetscMatrix>::IndexArray &permutation) {
        if (in.comm().size() == 1) {
            return false;
        }

        assert(in.rows() == in.cols());

        partitioning.resize(in.local_rows(), 0);

        if (!decompose(in, in.comm().size(), &partitioning[0])) {
            return false;
        }

        permutation.resize(in.local_rows(), 0);

        if (!partitions_to_permutations(in, &partitioning[0], &permutation[0])) {
            return false;
        }

        IS is = nullptr;
        PetscErrorCode err =
            ISCreateGeneral(in.comm().raw_comm(), in.local_rows(), &permutation[0], PETSC_USE_POINTER, &is);

        if (err != 0) {
            return false;
        }

        // Destroy because a new matrix is created below!
        out.destroy();
        err = MatPermute(in.raw_type(), is, is, &out.raw_type());

        ISDestroy(&is);
        return err == 0;
    }

}  // namespace utopia
