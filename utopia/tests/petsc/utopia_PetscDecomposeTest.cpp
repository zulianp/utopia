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
        UTOPIA_RUN_TEST(parmetis_partitions_to_permutations);
#endif
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

        std::stringstream ss;
        for (auto tag : decomposition) {
            ss << tag << "\n";
        }

        comm.synched_print(ss.str(), utopia::out().stream());
    }

    void parmetis_partitions_to_permutations() {
        auto &&comm = Comm::get_default();
        int mult = comm.rank() + 1;
        int n_local = mult * 3;

        Matrix mat;
        mat.sparse(layout(comm, n_local, n_local, Traits::determine(), Traits::determine()), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<SizeType> decomposition(n_local, -1);
        utopia_test_assert(parallel_decompose(mat, comm.size(), &decomposition[0]));

        IndexArray idx(n_local, -1);
        partitions_to_permutations(mat, &decomposition[0], &idx[0]);

        std::stringstream ss;
        for (auto i : idx) {
            ss << i << "\n";
        }

        comm.synched_print(ss.str(), utopia::out().stream());
    }
#endif
};

void petsc_decompose_test() {
#ifdef UTOPIA_WITH_PETSC
    DecomposeTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(petsc_decompose_test);
