#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_petsc_Decompose.hpp"

#include "utopia_petsc_SchurComplement.hpp"

#include "utopia_petsc_RebalancedSolver.hpp"

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

#ifdef UTOPIA_WITH_METIS
        UTOPIA_RUN_TEST(metis_decompose);
#endif

#ifdef UTOPIA_WITH_PARMETIS
        UTOPIA_RUN_TEST(parmetis_decompose);
        UTOPIA_RUN_TEST(parmetis_rebalance);
        UTOPIA_RUN_TEST(petsc_rebalanced_solver);
        UTOPIA_RUN_TEST(petsc_rebalanced_solver_file);
#endif
    }

#ifdef UTOPIA_WITH_METIS
    void metis_decompose() {
        Matrix mat;
        mat.sparse(serial_layout(20, 20), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<int> decomposition(mat.rows(), -1);
        utopia_test_assert(decompose(mat, 3, &decomposition[0]));

        //  for (auto tag : decomposition) {
        //            std::cout << tag << "\n";
        //}
    }
#endif

#ifdef UTOPIA_WITH_PARMETIS
    void parmetis_decompose() {
        auto &&comm = Comm::get_default();
        if (comm.size() == 1) return;

        int mult = comm.rank() + 1;

        Matrix mat;
        mat.sparse(layout(comm, mult * 3, mult * 3, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<int> decomposition(mat.local_rows(), -1);
        utopia_test_assert(parallel_decompose(mat, comm.size(), &decomposition[0]));

        // std::stringstream ss;
        // for (auto tag : decomposition) {
        //     ss << tag << "\n";
        // }

        // comm.synched_print(ss.str(), utopia::out().stream());
    }

    void parmetis_rebalance() {
        auto &&comm = Comm::get_default();
        if (comm.size() == 1) return;

        int mult = comm.rank() + 1;
        int n_local = mult * 3;

        Matrix mat;
        mat.sparse(layout(comm, n_local, n_local, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        Matrix rebalanced;
        std::vector<int> partitioning;
        IndexArray permutation;

        utopia_test_assert(rebalance(mat, rebalanced, partitioning, permutation));

        // disp(rebalanced);

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

    void petsc_rebalanced_solver() {
        auto &&comm = Comm::get_default();

        int n = 55;

        auto vl = layout(comm, Traits::decide(), n);
        auto ml = layout(comm, Traits::decide(), Traits::decide(), n, n);

        Matrix A;
        A.sparse(ml, 3, 2);
        Vector b(vl, 1.);
        Vector x(vl, 0.);

        assemble_laplacian_1D(A, true);

        {
            Range r = row_range(A);
            Write<Vector> w_b(b);

            if (r.inside(0)) {
                b.set(0, 0.);
            }

            if (r.inside(n - 1)) {
                b.set(n - 1, 0.);
            }
        }

        RebalancedSolver solver;

        auto p = utopia::param_list(utopia::param("inner_solver", utopia::param_list(utopia::param("verbose", true))));
        solver.read(p);

        solver.update(make_ref(A));
        solver.apply(b, x);

        int UTOPIA_DEBUG_RANDOM_DECOMPOSITION = 0;
        UTOPIA_READ_ENV(UTOPIA_DEBUG_RANDOM_DECOMPOSITION, atoi);

        // Test for multiple calls
        if (UTOPIA_DEBUG_RANDOM_DECOMPOSITION) {
            for (int i = 0; i < 10; i++) {
                x *= 0.999;
                solver.update(make_ref(A));
                solver.apply(b, x);
            }
        }
    }

    void petsc_rebalanced_solver_file() {
        auto &&comm = Comm::get_default();

        // Path path = Utopia::instance().get("data_path");
        // path = path / "forQR/A";

        Path path = "/Users/patrickzulian/Desktop/cloud/owncloud_HSLU/Patrick/2023/Cases/FP70/mat_raw/rowptr.raw";

        Matrix A;
        A.read(comm.get(), path);
        // A.read(comm.get(), "cippo/rowptr.raw");

        if (0) {
            // Check symmetry
            Matrix ApAt = A;
            ApAt.transform_values([](const Scalar) -> Scalar { return 1; });
            ApAt -= transpose(ApAt);
            double sum_ApAt = sum(abs(ApAt));

            printf("symmetric test %g == 0\n", sum_ApAt);
        }

        A.comm().synched_print(std::to_string(A.local_rows()));

        // A.write("AA.bin");

        auto vl = row_layout(A);
        Vector b(vl, 1.);
        Vector x(vl, 0.);

        // if (!A.comm().rank()) {
        //     utopia::out() << "dofs: " << A.rows() << "\n";
        // }

        // Temporary
        auto sp = utopia::param_list(utopia::param("type", "ksp"),
                                     utopia::param("ksp_type", "bcgstab"),
                                     utopia::param("pc_type", "bjacobi"),
                                     utopia::param("atol", "1e-8"),
                                     utopia::param("rtol", "1e-18"),
                                     utopia::param("stol", "1e-16"),
                                     utopia::param("verbose", true));

        auto p = utopia::param_list(
            // utopia::param("block_size", (comm.size() == 2 || comm.size() == 4)? 2 : 1),
            utopia::param("block_size", 4),
            // utopia::param("block_size", 1),
            utopia::param("inner_solver", std::move(sp)));

        RebalancedSolver solver;
        solver.read(p);

        // OmniLinearSolver<Matrix, Vector> solver;
        // p.get("inner_solver", solver);

        solver.update(make_ref(A));
        solver.apply(b, x);

        Vector r = b - A * x;
        Scalar r_norm = norm2(r);

        if (!r.comm().rank()) {
            utopia::out() << "r_norm: " << r_norm << "\n";
        }

        int UTOPIA_DEBUG_RANDOM_DECOMPOSITION = 0;
        UTOPIA_READ_ENV(UTOPIA_DEBUG_RANDOM_DECOMPOSITION, atoi);

        // Test for multiple calls
        if (UTOPIA_DEBUG_RANDOM_DECOMPOSITION) {
            for (int i = 0; i < 10; i++) {
                x *= 0.9999999;
                solver.update(make_ref(A));
                solver.apply(b, x);
            }
        }
    }
#endif

    void handcoded_partitions_to_permutations() {
        auto &&comm = Comm::get_default();

        if (comm.size() == 1) return;

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
        int n_local = 10;
        int n_eliminated = n_local - 2;

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
            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            // MPRGP<Matrix, Vector> cg;
            // cg.verbose(true);
            cg.apply_gradient_descent_step(true);

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
            ConjugateGradient<Matrix, Vector, HOMEMADE> og_cg;
            // MPRGP<Matrix, Vector> og_cg;
            og_cg.verbose(false);
            og_cg.apply_gradient_descent_step(true);

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
