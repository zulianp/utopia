#ifndef UTOPIA_PETSC_TYPES_HPP
#define UTOPIA_PETSC_TYPES_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    inline void convert(const Vec &petsc_vec, PetscVector &utopia_vec)
    {
        utopia_vec.convert_from(petsc_vec);
    }

    inline void convert(const PetscVector &utopia_vec, Vec &petsc_vec)
    {
        utopia_vec.convert_to(petsc_vec);
    }

    inline void convert(const Mat &petsc_mat, PetscMatrix &utopia_mat)
    {
        utopia_mat.convert_from(petsc_mat);
    }

    inline void convert(const PetscMatrix &utopia_mat, Mat &petsc_mat)
    {
        utopia_mat.convert_to(petsc_mat);
    }

    inline void wrap(Mat &m, PetscMatrix &utopia_mat) {
        utopia_mat.wrap(m);
    }

    inline void wrap(Vec &vec, PetscVector &utopia_vec)
    {
        utopia_vec.wrap(vec);
    }

    inline int comm_size(const PetscVector &t)
    {
        auto comm = t.comm().get();
        int ret;
        MPI_Comm_size(comm, &ret);
        return ret;
    }

    inline int comm_rank(const PetscVector &t)
    {
        auto comm = t.comm().get();
        int ret;
        MPI_Comm_rank(comm, &ret);
        return ret;
    }

    inline int comm_size(const PetscMatrix &t)
    {
        auto comm = t.comm().get();
        int ret;
        MPI_Comm_size(comm, &ret);
        return ret;
    }

    inline int comm_rank(const PetscMatrix &t)
    {
        auto comm = t.comm().get();
        int ret;
        MPI_Comm_rank(comm, &ret);
        return ret;
    }

    inline void synchronize(PetscVector &t)
    {
        t.update_ghosts();
    }

    ///////////////////////////////////////////////////////////////////////////////

}

#endif //UTOPIA_UTOPIA_PETSC_TYPES_HPP

