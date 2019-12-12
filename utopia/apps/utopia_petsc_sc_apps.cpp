
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
//include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Core.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_petsc_Matrix.hpp"

#include <cmath>

namespace utopia {

    static void petsc_dm_app()
    {
        using SizeType = Traits<PetscVector>::SizeType;

        PetscCommunicator world;
        PetscDM dm(
            world,
            {3, 3},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        dm.describe();

        PetscMatrix mat;
        dm.create_matrix(mat);

        disp(mat);

        dm.each_element([](const SizeType &i, const PetscDM::Elem &e) {

        });
    }

    UTOPIA_REGISTER_APP(petsc_dm_app);
}

#endif //WITH_TRILINOS
