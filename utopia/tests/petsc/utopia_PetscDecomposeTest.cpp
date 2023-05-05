#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_petsc_Decompose.hpp"

#include "utopia_petsc_SchurComplement.hpp"

#include "petscmat.h"

using namespace utopia;

template <class Matrix, class Vector>
class DecomposeTest {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;
    using IndexArray = typename Traits::IndexArray;
    using Comm = typename Traits::Communicator;

    void run() {
#ifdef UTOPIA_WITH_METIS
        UTOPIA_RUN_TEST(metis_decompose);
#endif

#ifdef UTOPIA_WITH_PARMETIS
        UTOPIA_RUN_TEST(parmetis_decompose);
        UTOPIA_RUN_TEST(parmetis_rebalance);
#endif

        UTOPIA_RUN_TEST(handcoded_partitions_to_permutations);
// FIXME(zulianp) this test causes Segmentation Violation when utopia is built with gpu support
        if constexpr (Traits::Backend == PETSC) {
#ifdef PETSC_HAVE_CUDA
            utopia_warning("Skipping schur_complement");
#else
            UTOPIA_RUN_TEST(schur_complement);
#endif
        } else {
            UTOPIA_RUN_TEST(schur_complement);
        }
    }

#ifdef UTOPIA_WITH_METIS
    void metis_decompose() {
        Matrix mat;
        mat.sparse(serial_layout(20, 20), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<SizeType> decomposition(mat.rows(), -1);
        utopia_test_assert(decompose(mat, 3, &decomposition[0]));

        //  for (auto tag : decomposition) {
        //            std::cout << tag << "\n";
        //}
    }
#endif

#ifdef UTOPIA_WITH_PARMETIS
    void parmetis_decompose() {
        auto &&comm = Comm::get_default();
        int mult = comm.rank() + 1;

        Matrix mat;
        mat.sparse(layout(comm, mult * 3, mult * 3, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<SizeType> decomposition(mat.local_rows(), -1);
        utopia_test_assert(parallel_decompose(mat, comm.size(), &decomposition[0]));

        // std::stringstream ss;
        // for (auto tag : decomposition) {
        //     ss << tag << "\n";
        // }

        // comm.synched_print(ss.str(), utopia::out().stream());
    }

    void parmetis_rebalance() {
        auto &&comm = Comm::get_default();
        int mult = comm.rank() + 1;
        int n_local = mult * 3;

        Matrix mat;
        mat.sparse(layout(comm, n_local, n_local, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        Matrix rebalanced;
        IndexArray partitioning, permutation;

        utopia_test_assert(rebalance(mat, rebalanced, partitioning, permutation));

#if 0
        {
            auto rr = mat.row_ranges();

            if (comm.rank() == 0) {
                for (auto r : rr) {
                    std::cout << r << " ";
                }

                std::cout << "\n";
            }
        }

        {
            auto rr = rebalanced.row_ranges();

            if (comm.rank() == 0) {
                for (auto r : rr) {
                    std::cout << r << " ";
                }

                std::cout << "\n";
            }
        }

        disp(rebalanced);
#endif
    }
#endif

    void handcoded_partitions_to_permutations() {
        auto &&comm = Comm::get_default();
        int mult = comm.rank() + 1;
        int n_local = mult * 2;

        Matrix mat;
        mat.sparse(layout(comm, n_local, n_local, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<int> decomposition(n_local, 0);

        for (int i = 0; i < n_local; ++i) {
            // shift entries to 3 ranks above in a circular way
            decomposition[i] = (comm.rank() + 2) % comm.size();
        }

        IndexArray partitioning;
        utopia_test_assert(partitions_to_permutations(mat, &decomposition[0], partitioning));

        Matrix redistributed;
        redistribute_from_permutation(mat, partitioning, redistributed);

        std::vector<int> inverse_decomposition;
        inverse_partition_mapping(comm.size(), mat.row_ranges(), partitioning, inverse_decomposition);

        IndexArray inverse_partitioning;
        utopia_test_assert(partitions_to_permutations(redistributed, &inverse_decomposition[0], inverse_partitioning));

#if 0
        std::stringstream ss;
        ss << "inverse_decomposition: ";
        for (auto sd : inverse_decomposition) {
            ss << sd << " ";
        }

        ss << "\n";

        ss << "partitioning: ";
        for (auto sd : partitioning) {
            ss << sd << " ";
        }

        ss << "\n";

        ss << "inverse_partitioning: ";
        for (auto sd : inverse_partitioning) {
            ss << sd << " ";
        }

        ss << "\n";

        comm.synched_print(ss.str(), utopia::out().stream());
#endif

        Matrix inverse_redistributed;
        redistribute_from_permutation(redistributed, inverse_partitioning, inverse_redistributed);

        // disp(mat);

        // disp(redistributed);

        // disp(inverse_redistributed);

        utopia_test_assert(approxeq(inverse_redistributed, mat));
    }

    void schur_complement() {
        auto &&comm = Comm::get_default();
        int n_local = 1000;
        int n_eliminated = n_local - 1;

        Matrix mat;
        mat.sparse(layout(comm, n_local, n_local, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat, true);

        auto rr = mat.row_range();
        Vector rhs(row_layout(mat), 1);

        int eliminate_offset = 1;
        if (rr.inside(0)) {
            Write<Vector> w(rhs);
            rhs.set(0, 0);
        }

        if (rr.inside(mat.rows() - 1)) {
            Write<Vector> w(rhs);
            rhs.set(mat.rows() - 1, 0);
        }

        IndexArray idx(n_eliminated);
        for (int i = 0; i < n_eliminated; ++i) {
            idx[i] = rr.begin() + eliminate_offset + i;
        }

        Chrono sc_time;
        sc_time.start();

        SchurComplement<Matrix> sc;
        sc.initialize_from_selection(mat, idx);

        Vector x(layout(rhs), 0);
        Vector rhs_eliminated, rhs_sc;
        sc.apply_righthand_side(rhs, rhs_eliminated, rhs_sc);
        {
            // ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            MPRGP<Matrix, Vector> cg;
            cg.verbose(true);
            // cg.apply_gradient_descent_step(true);

            cg.max_it(40000);
            cg.atol(1e-50);
            cg.stol(1e-18);
            cg.rtol(1e-18);

            Vector x_sc(layout(rhs_sc), 0);
            cg.solve(sc, rhs_sc, x_sc);
            sc.finalize_from_eliminated_rhs(rhs_eliminated, x_sc, x);
        }

        sc_time.stop();

        Scalar norm_sc = norm2(rhs - mat * x);

        Chrono og_time;
        og_time.start();

        // Oracle
        Vector og_x(row_layout(mat), 0);
        {
            // ConjugateGradient<Matrix, Vector, HOMEMADE> og_cg;
            MPRGP<Matrix, Vector> og_cg;
            og_cg.verbose(true);
            // og_cg.apply_gradient_descent_step(true);

            og_cg.max_it(40000);
            og_cg.atol(1e-8);
            og_cg.stol(1e-18);
            og_cg.rtol(1e-18);

            og_cg.solve(mat, rhs, og_x);
        }

        og_time.stop();

        Scalar norm_og = norm2(rhs - mat * og_x);

        std::stringstream ss;
        ss << "\nSchur " << norm_sc << "\n" << sc_time;
        ss << "OG " << norm_og << "\n" << og_time;

        mat.comm().root_print(ss.str());

        // Vector diff = og_x - x;
        // Scalar diff_norm = norm2(diff);

        // disp(diff_norm);
        // utopia_test_assert(diff_norm < 1e-6);
    }
};

void petsc_decompose_test() {
#ifdef UTOPIA_WITH_PETSC
    DecomposeTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(petsc_decompose_test);
