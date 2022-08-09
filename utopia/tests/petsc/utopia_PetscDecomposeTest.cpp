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
    using Comm = typename Traits::Communicator;

    void run() {
#ifdef UTOPIA_WITH_METIS
        UTOPIA_RUN_TEST(metis_decompose);
#endif
    }

    void metis_decompose() {
#ifdef UTOPIA_WITH_METIS
        Matrix mat;
        mat.sparse(serial_layout(20, 20), 3, 3);
        assemble_laplacian_1D(mat);

        std::vector<SizeType> decomposition(mat.rows(), -1);
        utopia_test_assert(decompose(mat, 3, &decomposition[0]));

        for (auto tag : decomposition) {
            std::cout << tag << "\n";
        }

#endif
    }
};

void petsc_decompose_test() {
#ifdef UTOPIA_WITH_PETSC
    DecomposeTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(petsc_decompose_test);
