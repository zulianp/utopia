
#include <cstring>
#include <string>
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_Operator.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {

    PetscErrorCode UtopiaOperatorApplyShell(Mat m, Vec X, Vec Y) {
        Operator<PetscVector> *ut_operator;
        MatShellGetContext(m, &ut_operator);

        PetscVector x_u, y_u;
        // convert(X, x_u);
        // convert(Y, y_u);
        wrap(X, x_u);
        wrap(Y, y_u);

        ut_operator->apply(x_u, y_u);
        // convert(y_u, Y);

        return (0);
    }

}  // namespace utopia
