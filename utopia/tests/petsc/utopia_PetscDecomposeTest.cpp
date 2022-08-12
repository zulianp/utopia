#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_petsc_Decompose.hpp"

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
    }
#endif

    void handcoded_partitions_to_permutations() {
        auto &&comm = Comm::get_default();
        int mult = comm.rank() + 1;
        int n_local = mult * 2;

        Matrix mat;
        mat.sparse(layout(comm, n_local, n_local, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<SizeType> decomposition(n_local, 0);

        auto rr = mat.row_range();
        SizeType n_rows = mat.rows();
        for (int i = 0; i < n_local; ++i) {
            decomposition[i] = (comm.rank() + 1) % comm.size();
        }

        IndexArray partitioning;
        utopia_test_assert(partitions_to_permutations(mat, &decomposition[0], partitioning));

        Matrix redistributed;
        redistribute_from_permutation(mat, partitioning, redistributed);

        IndexArray inverse_decomposition;
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
};

void petsc_decompose_test() {
#ifdef UTOPIA_WITH_PETSC
    DecomposeTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(petsc_decompose_test);
