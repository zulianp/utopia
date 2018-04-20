
#ifndef UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP
#define UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

namespace utopia 
{
    class PetscMatrix;
    class PetscSparseMatrix;
    class PetscSerialSparseMatrix;

    class PetscVector;
    class PetscSerialVector;

    template<class Matrix, class Vector>
    class TaoSolver;
}

#endif //UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

