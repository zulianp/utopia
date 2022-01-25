#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "petscmat.h"

using namespace utopia;

template <class Matrix, class Vector>
class CommStructTest {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;
    using Comm = typename Traits::Communicator;

    void run() { UTOPIA_RUN_TEST(test_comm_struct); }

    void test_comm_struct() {
        auto &&comm = Comm::get_default();
        Matrix mat;
        mat.sparse(square_matrix_layout(layout(comm, 4, 4 * comm.size())), 3, 3);
        assemble_laplacian_1D(mat);

        disp(mat);

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 15, 3)

        Vec lvec;
        PetscInt *colmap;
        VecScatter scatt;

        MatGetCommunicationStructs(mat.raw_type(), &lvec, &colmap, &scatt);
        // disp(colmap);

        VecScatterView(scatt, PETSC_VIEWER_STDOUT_WORLD);
#endif  // UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN

        // MatGetGhosts
        PetscInt nghosts;
        const PetscInt *ghosts;
        MatGetGhosts(mat.raw_type(), &nghosts, &ghosts);

        std::stringstream ss;
        ss << "nghosts: " << nghosts << "\n";
        for (SizeType i = 0; i < nghosts; ++i) {
            ss << ghosts[i] << " ";
        }

        mat.comm().synched_print(ss.str());

        // PetscErrorCode MatGetOwnershipRanges(Mat mat,const PetscInt **ranges)

        // VecScatter scatter;
        // MatScatterGetVecScatter(mat.raw_type(), &scatter);
        // VecScatterView(scatter, PETSC_VIEWER_STDOUT_WORLD);
    }
};

void petsc_comm_struct() {
#ifdef UTOPIA_WITH_PETSC
    CommStructTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(petsc_comm_struct);
