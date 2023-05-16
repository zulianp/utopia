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
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_Utils.hpp"

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

#ifndef UTOPIA_WITH_PARMETIS

    bool parallel_decompose(const PetscMatrix &, const int, int *) { return false; }

#else

    bool parallel_decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions) {
        auto rrs = matrix.row_ranges();

        int comm_size = matrix.comm().size();

        idx_t *vtxdist = nullptr;
        if (sizeof(idx_t) != sizeof(PetscInt)) {
            vtxdist = (idx_t *)malloc((comm_size + 1) * sizeof(idx_t));
            for (int r = 0; r <= comm_size; r++) {
                vtxdist[r] = rrs[r];
            }
        } else {
            vtxdist = (idx_t *)&rrs[0];
        }

        const idx_t row_offset = vtxdist[matrix.comm().rank()];

        idx_t ncon = 1;
        idx_t *vwgt = nullptr;
        // idx_t* vsize = nullptr;
        idx_t *adjwgt = nullptr;
        idx_t wgtflag = 0;
        idx_t numflag = 0;
        idx_t nparts = num_partitions;

        real_t ubvec[1] = {1.05};
        idx_t options[3] = {0};
        // idx_t objval = -1;
        idx_t edgecut;
        idx_t *parts = partitions;
        MPI_Comm comm = matrix.comm().raw_comm();

        Mat d, o;

        const PetscInt *colmap;
        MatMPIAIJGetSeqAIJ(matrix.raw_type(), &d, &o, &colmap);

        PetscCrsView d_view(d);
        PetscCrsView o_view(o);

        std::vector<idx_t> xadj(d_view.rows() + 1, 0);
        std::vector<idx_t> adjncy(d_view.colidx().size() + o_view.colidx().size(), -1);
        std::vector<real_t> tpwgts(num_partitions, 1. / num_partitions);
        std::vector<idx_t> actual_vwgts(d_view.rows(), 0);

        auto cr = matrix.col_range();

        PetscInt local_rows = d_view.rows();
        PetscInt col_offset = cr.begin();

        auto d_xadj = d_view.row_ptr();
        xadj[0] = 0;

        for (PetscInt i = 0; i < local_rows; ++i) {
            auto d_row = d_view.row(i);
            auto o_row = o_view.row(i);

            xadj[i + 1] = xadj[i];
            for (PetscInt k = 0; k < d_row.length; ++k) {
                const idx_t col = d_row.colidx(k);
                if (col_offset + col != row_offset + i) {
                    adjncy[xadj[i + 1]++] = col_offset + col;
                }
            }

            for (PetscInt k = 0; k < o_row.length; ++k) {
                const idx_t col = colmap[o_row.colidx(k)];
                if (col != row_offset + i) {
                    adjncy[xadj[i + 1]++] = col;
                }
            }
        }

        for (PetscInt i = 0; i < local_rows; ++i) {
            actual_vwgts[i] = xadj[i + 1] - xadj[i];
        }

        // vwgt = &actual_vwgts[0];
        // wgtflag = 2;

        vwgt = 0;
        wgtflag = 0;

        int ret = ParMETIS_V3_PartKway(vtxdist,     // 0
                                       &xadj[0],    // 1
                                       &adjncy[0],  // 2
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

        if (sizeof(idx_t) != sizeof(PetscInt)) {
            free(vtxdist);
        }

        if (ret == METIS_OK) {
            return true;
        } else {
            return false;
        }
    }
#endif

    bool partitions_to_permutations(const Communicator &comm,
                                    const ArrayView<const PetscInt> &rrs,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index,
                                    std::vector<int> &rpartitions,
                                    Traits<PetscMatrix>::IndexArray &rindex) {
        const int comm_size = comm.size();

        PetscInt r_begin = rrs[comm.rank()];

        std::vector<PetscInt> send_row_count(comm_size, 0);
        std::vector<PetscInt> recv_row_count(comm_size, 0);

        const int local_rows = rrs[comm.rank() + 1] - r_begin;

        std::vector<PetscInt> sendbuf(local_rows, 0);
        std::vector<int> sdispls(comm_size + 1, 0);
        std::vector<int> sendcounts(comm_size, 0);

        rindex.resize(local_rows);

        {  // Fill send buffers!

            for (int i = 0; i < local_rows; ++i) {
                int part = partitions[i];
                assert(part < comm.size());
                assert(part >= 0);

                rindex[i] = sdispls[part + 1];
                ++sdispls[part + 1];
            }

            for (int i = 0; i < comm_size; ++i) {
                sdispls[i + 1] += sdispls[i];
            }

            for (int i = 0; i < local_rows; ++i) {
                int part = partitions[i];
                int idx = sdispls[part] + sendcounts[part]++;
                sendbuf[idx] = r_begin + i;
            }
        }

        std::vector<int> rdispls(comm_size + 1, 0);
        std::vector<int> recvcounts(comm_size, 0);
        std::vector<int> roffsets(comm_size, 0);
        std::vector<PetscInt> sizes(comm_size, 0);
        PetscInt incoming = 0;

        {  // Fill recv buffers
            // exchange send/receive permuatation info
            MPI_Alltoall(sendcounts.data(),
                         1,
                         utopia::MPIType<int>::value(),
                         recvcounts.data(),
                         1,
                         utopia::MPIType<int>::value(),
                         comm.raw_comm());

            for (int r = 0; r < comm_size; ++r) {
                // Rows that are kept local are also counted!
                incoming += recvcounts[r];
                rdispls[r + 1] = rdispls[r] + recvcounts[r];
            }

            MPI_Alltoall(rdispls.data(),
                         1,
                         utopia::MPIType<int>::value(),
                         roffsets.data(),
                         1,
                         utopia::MPIType<int>::value(),
                         comm.raw_comm());

            MPI_Allgather(
                &incoming, 1, MPIType<PetscInt>::value(), sizes.data(), 1, MPIType<PetscInt>::value(), comm.raw_comm());

            PetscInt proc_offset = 0;
            for (int r = 0; r < comm_size; r++) {
                roffsets[r] += proc_offset;
                proc_offset += sizes[r];
            }
        }

        index.resize(incoming, 0);

        MPI_Alltoallv(sendbuf.data(),
                      sendcounts.data(),
                      sdispls.data(),
                      MPIType<PetscInt>::value(),
                      index.data(),
                      recvcounts.data(),
                      rdispls.data(),
                      MPIType<PetscInt>::value(),
                      comm.raw_comm());

        rpartitions.resize(incoming, -1);

        for (SizeType r = 0; r < comm_size; r++) {
            SizeType offset = rdispls[r];
            SizeType count = recvcounts[r];

            for (SizeType i = 0; i < count; i++) {
                rpartitions[offset + i] = r;
            }
        }

        for (int i = 0; i < local_rows; ++i) {
            int part = partitions[i];
            assert(part < comm.size());
            assert(part >= 0);

            rindex[i] += roffsets[part];
        }

        return true;
    }

    bool partitions_to_permutations(const Communicator &comm,
                                    const ArrayView<const PetscInt> &rrs,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index) {
        const int comm_size = comm.size();

        PetscInt r_begin = rrs[comm.rank()];

        std::vector<PetscInt> send_row_count(comm_size, 0);
        std::vector<PetscInt> recv_row_count(comm_size, 0);

        const int local_rows = rrs[comm.rank() + 1] - r_begin;

        std::vector<PetscInt> sendbuf(local_rows, 0);
        std::vector<int> sdispls(comm_size + 1, 0);
        std::vector<int> sendcounts(comm_size, 0);

        {  // Fill send buffers!

            for (int i = 0; i < local_rows; ++i) {
                int part = partitions[i];
                assert(part < comm.size());
                assert(part >= 0);
                ++sdispls[part + 1];
            }

            for (int i = 0; i < comm_size; ++i) {
                sdispls[i + 1] += sdispls[i];
            }

            for (int i = 0; i < local_rows; ++i) {
                int part = partitions[i];
                int idx = sdispls[part] + sendcounts[part]++;
                sendbuf[idx] = r_begin + i;
            }
        }

        std::vector<int> rdispls(comm_size + 1, 0);
        std::vector<int> recvcounts(comm_size, 0);
        PetscInt incoming = 0;

        {  // Fill recv buffers
            // exchange send/receive permuatation info
            MPI_Alltoall(sendcounts.data(),
                         1,
                         utopia::MPIType<int>::value(),
                         recvcounts.data(),
                         1,
                         utopia::MPIType<int>::value(),
                         comm.raw_comm());

            for (int r = 0; r < comm_size; ++r) {
                // Rows that are kept local are also counted!
                incoming += recvcounts[r];
                rdispls[r + 1] = rdispls[r] + recvcounts[r];
            }
        }

        index.resize(incoming, 0);

        MPI_Alltoallv(sendbuf.data(),
                      sendcounts.data(),
                      sdispls.data(),
                      MPIType<PetscInt>::value(),
                      index.data(),
                      recvcounts.data(),
                      rdispls.data(),
                      MPIType<PetscInt>::value(),
                      comm.raw_comm());

#if 0
            {
                std::stringstream ss;

                ss << "partitions: ";

                for (int i = 0; i < local_rows; ++i) {
                    ss << partitions[i] << " ";
                }

                ss << "\n";

                ss << "sdispls: ";
                for (auto sd : sdispls) {
                    ss << sd << " ";
                }

                ss << "\n";

                ss << "sendcounts: ";
                for (auto sc : sendcounts) {
                    ss << sc << " ";
                }

                ss << "\n";

                ss << "recvcounts: ";
                for (auto rc : recvcounts) {
                    ss << rc << " ";
                }

                ss << "\n";

                ss << "rdispls: ";
                for (auto rd : rdispls) {
                    ss << rd << " ";
                }

                ss << "\n";

                ss << "sendbuf: ";
                for (auto idx : sendbuf) {
                    ss << idx << " ";
                }

                ss << "\n";

                ss << "index: ";
                for (auto idx : index) {
                    ss << idx << " ";
                }

                ss << "\n";

                comm.synched_print(ss.str(), utopia::out().stream());
            }
#endif

        return true;
    }

    bool partitions_to_permutations(const PetscMatrix &matrix,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index) {
        partitions_to_permutations(matrix.comm(), matrix.row_ranges(), partitions, index);
        return true;
    }

    bool inverse_partition_mapping(const int comm_size,
                                   const ArrayView<const PetscInt> &original_row_ranges,
                                   const Traits<PetscMatrix>::IndexArray &permutation,
                                   std::vector<int> &partitioning) {
        partitioning.resize(permutation.size());
        int n_local = original_row_ranges[1] - original_row_ranges[0];

        PetscInt n_indices = permutation.size();
        for (PetscInt i = 0; i < n_indices; ++i) {
            partitioning[i] = find_rank(comm_size, n_local, &original_row_ranges[0], permutation[i]);

            // printf("%d -> %d\n", permutation[i], partitioning[i]);

            assert(partitioning[i] < comm_size);
            assert(partitioning[i] >= 0);
        }

        return true;
    }

    bool redistribute_from_permutation(const PetscMatrix &in,
                                       const Traits<PetscMatrix>::IndexArray &permutation,
                                       PetscMatrix &out,
                                       MatReuse reuse) {
        IS is = nullptr;
        PetscErrorCode err =
            ISCreateGeneral(in.comm().raw_comm(), permutation.size(), &permutation[0], PETSC_USE_POINTER, &is);

        if (err != 0) {
            return false;
        }

        out.destroy();  // Destroy because a new matrix is created below!
        err = MatCreateSubMatrix(in.raw_type(), is, is, reuse, &out.raw_type());

        ISDestroy(&is);
        return true;
    }

    // bool redistribute_columns_from_permutation(const PetscMatrix &in,
    //                                            const Traits<PetscMatrix>::IndexArray &permutation,
    //                                            PetscMatrix &out,
    //                                            MatReuse reuse) {
    //     IS is_col = nullptr;
    //     PetscErrorCode err =
    //         ISCreateGeneral(in.comm().raw_comm(), permutation.size(), &permutation[0], PETSC_USE_POINTER, &is_col);

    //     if (err != 0) {
    //         return false;
    //     }

    //     IS is_row = nullptr;

    //     out.destroy();  // Destroy because a new matrix is created below!
    //     err = MatCreateSubMatrix(in.raw_type(), is_col, is_col, reuse, &out.raw_type());

    //     ISDestroy(&is_col);
    //     return true;
    // }

    bool redistribute_from_permutation(const PetscVector &in,
                                       const Traits<PetscMatrix>::IndexArray &permutation,
                                       PetscVector &out) {
        in.select(permutation, out);
        return true;
    }

    bool rebalance(const PetscMatrix &in,
                   PetscMatrix &out,
                   std::vector<int> &partitioning,
                   Traits<PetscMatrix>::IndexArray &permutation) {
        if (in.comm().size() == 1) {
            return false;
        }

        assert(in.rows() == in.cols());

        partitioning.resize(in.local_rows(), 0);

        if (!parallel_decompose(in, in.comm().size(), &partitioning[0])) {
            return false;
        }

        if (!partitions_to_permutations(in, &partitioning[0], permutation)) {
            return false;
        }

        return redistribute_from_permutation(in, permutation, out);
    }

    bool rebalance(const PetscMatrix &in,
                   PetscMatrix &out,
                   std::vector<int> &partitioning,
                   Traits<PetscMatrix>::IndexArray &permutation,
                   std::vector<int> &r_partitioning,
                   Traits<PetscMatrix>::IndexArray &r_permutation) {
        if (in.comm().size() == 1) {
            return false;
        }

        assert(in.rows() == in.cols());

        partitioning.resize(in.local_rows(), 0);

        if (!parallel_decompose(in, in.comm().size(), &partitioning[0])) {
            return false;
        }

        if (!partitions_to_permutations(
                in.comm(), in.row_ranges(), &partitioning[0], permutation, r_partitioning, r_permutation)) {
            return false;
        }

        bool ok = redistribute_from_permutation(in, permutation, out);
        return ok;
    }

}  // namespace utopia
