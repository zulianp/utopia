#include "utopia_petsc_Eval_Parallel.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_quirks.hpp"

namespace utopia {

    void build_local_redistribute(
        const PetscVector &x_from,
        const PetscVector &shape_vec,
        PetscVector &result)
    {
        if(x_from.comm().size() == 1)
        {
            result = x_from;
        } else {
            PetscVector x_to = shape_vec;
            IS is;
            VecScatter newctx;

            PetscInt ed, st;
            VecGetOwnershipRange(x_to.raw_type() ,&st, &ed);
            ISCreateStride(PETSC_COMM_SELF, (ed - st), st, 1, &is);

            UtopiaVecScatterCreate(x_from.raw_type(), is, x_to.raw_type(), is, &newctx);

            VecScatterBegin(newctx, x_from.raw_type(), x_to.raw_type(), INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(newctx, x_from.raw_type(), x_to.raw_type(), INSERT_VALUES, SCATTER_FORWARD);


            result = x_to;

            VecScatterDestroy(&newctx);
            ISDestroy(&is);
        }
    }

}
