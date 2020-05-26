#ifndef UTOPIA_PETSC_LAYOUT_HPP
#define UTOPIA_PETSC_LAYOUT_HPP

#include "petscmat.h"
#include "petscvec.h"
#include "utopia_petsc_Traits.hpp"

namespace utopia {

    inline PetscTraits::Layout layout(Vec v) {
        PetscInt n_local, n_global;
        VecGetLocalSize(v, &n_local);
        VecGetSize(v, &n_global);

        MPI_Comm comm = PetscObjectComm((PetscObject)v);
        assert(comm != MPI_COMM_NULL);
        return layout(PetscTraits::Communicator(comm), n_local, n_global);
    }

    inline PetscTraits::MatrixLayout layout(Mat m) {
        PetscInt rows, cols;
        MatGetSize(m, &rows, &cols);

        PetscInt local_rows, local_cols;
        MatGetLocalSize(m, &local_rows, &local_cols);

        MPI_Comm comm = PetscObjectComm((PetscObject)m);
        assert(comm != MPI_COMM_NULL);
        return layout(PetscTraits::Communicator(comm), local_rows, local_cols, rows, cols);
    }

}  // namespace utopia

#endif  // UTOPIA_PETSC_LAYOUT_HPP
