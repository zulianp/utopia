#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include <petscdm.h>

namespace utopia {
    class PetscDMBase {
    public:
        virtual ~PetscDMBase() {}
        virtual void init(DM dm) = 0;

        template<class SizeType>
        static void dof_ownership_range(const DM &dm, SizeType &begin, SizeType &end)
        {
            Vec v;
            DMGetGlobalVector(dm, &v);
            VecGetOwnershipRange(v, &begin, &end);
            DMRestoreGlobalVector(dm, &v);
        }

        static PetscInt get_dimension(DM dm)
        {
            PetscInt ret;
            DMGetDimension(dm, &ret);
            return ret;
        }
    };
}

#endif //UTOPIA_PETSC_DM_HPP
