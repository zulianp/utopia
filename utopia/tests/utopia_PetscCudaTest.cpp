#include "utopia.hpp"
#include "utopia_PetscCudaTest.hpp"

namespace utopia {

#ifdef WITH_PETSC
    void petsc_cuda_init()
    {

    }

#endif //WITH_PETSC;

    void run_petsc_cuda_test() {
#ifdef WITH_PETSC
        UTOPIA_UNIT_TEST_BEGIN("PetscCudaTest");
        UTOPIA_RUN_TEST(petsc_cuda_init);
        UTOPIA_UNIT_TEST_END("PetscCudaTest");
#endif // WITH_PETSC
    }
}
