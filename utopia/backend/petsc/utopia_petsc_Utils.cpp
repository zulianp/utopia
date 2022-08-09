#include "utopia_petsc_Utils.hpp"

#include "utopia_petsc_Matrix_impl.hpp"

#include "utopia_Diag.hpp"
#include "utopia_Writable.hpp"

#include "utopia_Core.hpp"

#include "utopia_DeviceView.hpp"
#include "utopia_RowView.hpp"
#include "utopia_petsc_RowView.hpp"

namespace utopia {

    void optimize_nnz(PetscMatrix &A) {
        UTOPIA_TRACE_REGION_BEGIN("optimize_nnz");

        auto rr = row_range(A);
        auto cr = col_range(A);
        auto ls = local_size(A);
        auto gs = size(A);

        std::vector<PetscInt> d_nnz(rr.extent(), 0), o_nnz(rr.extent(), 0);
        A.read([&](const utopia::SizeType i, const utopia::SizeType j, const PetscScalar val) {
            if (std::abs(val) > 1e-18) {
                if (cr.inside(j)) {
                    ++d_nnz[i - rr.begin()];
                } else {
                    ++o_nnz[i - rr.begin()];
                }
            }
        });

        PetscMatrix A_opt;

        A_opt.matij_init(A.communicator(), A.type_override(), ls.get(0), ls.get(1), gs.get(0), gs.get(1), d_nnz, o_nnz);

        {
            Write<PetscMatrix> w_A(A_opt);
            A.read([&](const SizeType i, const SizeType j, const PetscScalar val) {
                if (std::abs(val) > 1e-18 || i == j) {
                    A_opt.set(i, j, val);
                }
            });
        }

        A = std::move(A_opt);

        UTOPIA_TRACE_REGION_END("optimize_nnz");
    }

    bool is_diagonally_dominant(const PetscMatrix &A) {
        PetscVector d = diag(A);
        PetscVector o(layout(d));

        {
            Write<PetscVector> w_o(o);
            A.read([&o](const SizeType i, const SizeType j, const PetscScalar val) {
                if (i != j) {
                    o.add(i, std::abs(val));
                }
            });
        }

        PetscVector diff = d - o;
        PetscScalar m = min(diff);
        return m > 0.;
    }

    void local_block_view(const PetscMatrix &mat, PetscMatrix &block) {
        Mat M;
        auto ierr = MatGetDiagonalBlock(mat.raw_type(), &M);
        assert(ierr == 0);
        UTOPIA_UNUSED(ierr);

        block.wrap(M);
        block.update_mirror();
    }

    UTOPIA_FUNCTION int find_rank(int comm_size, PetscInt n_local, const PetscInt *ranges, const Size_t global_id) {
        int rank = device::min(int(global_id * (float(comm_size) / n_local)), comm_size - 1);

        bool found = (global_id >= ranges[rank]) && (global_id < ranges[rank + 1]);

        while (!found) {
            if (global_id < ranges[rank]) {
                --rank;
            } else if (global_id >= ranges[rank + 1]) {
                ++rank;
            } else {
                assert(false);
            }

            assert(global_id >= 0);
            assert(global_id < ranges[comm_size]);
            assert(rank < comm_size);
            assert(rank >= 0);

            found = (global_id >= ranges[rank]) && (global_id < ranges[rank + 1]);
        }

        assert(found);
        return rank;
    }

    bool transpose_distro_is_strongly_imbalanced(const PetscMatrix &mat, double tol) {
        PetscVector w;
        compute_column_nnz_weight(mat, w);

        double w_min = min(w);
        double w_max = max(w);

        double diff = w_max - w_min;
        return diff > tol;
    }

    void compute_column_nnz_weight(const PetscMatrix &mat, PetscVector &weights) {
        using TraitsT = utopia::Traits<PetscMatrix>;

        using Comm = typename TraitsT::Communicator;
        using Scalar = typename TraitsT::Scalar;
        using SizeType = typename TraitsT::SizeType;

        auto &&comm = mat.comm();

        int n_local = mat.local_cols();
        int comm_size = comm.size();
        int comm_rank = comm.rank();

        weights.zeros(serial_layout(comm_size));

        {
            auto v_view = local_view_device(weights);

            auto cr = mat.col_ranges();
            mat.read(UTOPIA_LAMBDA(const SizeType, const SizeType j, const Scalar) {
                int r = find_rank(comm_size, n_local, &cr[0], j);
                v_view.add(r, 1);
            });

            comm.sum(comm_size, &v_view.array()[0]);
        }

        // if (comm.rank() == 0) {
        //     disp(weights);
        // }
    }

    // void compute_col_rebalancing(const PetscCommunicator &comm,
    //                              const PetscVector &weights,
    //                              const PetscInt *original_colranges,
    //                              PetscInt *colranges) {
    //     int comm_size = comm.size();
    //     auto w_view = local_view_device(weights);
    //     double total_work = sum(weights);
    //     double avg_work = total_work / comm_size;
    //     PetscInt n_cols = original_colranges[comm_size];

    //     for (int r = 0; r < comm_size; ++r) {
    //         double work = w_view.get(r);
    //         int n_local = original_colranges[r + 1] - original_colranges[r];
    //         double ratio = avg_work / work;
    //         int new_nlocal = std::ceil(n_local * ratio);
    //         colranges[r + 1] = std::min(n_cols, colranges[r] + new_nlocal);
    //     }
    // }

}  // namespace utopia
