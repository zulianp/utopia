
// #ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
// #define UTOPIA_PETSC_SPARSE_MATRIX_HPP

// #include "utopia_petsc_Matrix.hpp"
// #include <map>
// #include "petscmat.h"

// namespace utopia{

//     class PetscSparseMatrix : public PetscMatrix {
//     public:
//         virtual MatType type_override() const override
//         {
//             return MATAIJ;
//         }

//         virtual ~PetscSparseMatrix() {}
//     };

//     class PetscCuSparseMatrix : public PetscMatrix {
//     public:
//         virtual MatType type_override() const override
//         {
//             return MATAIJCUSPARSE;
//         }

//         virtual ~PetscCuSparseMatrix() {}
//     };

// }

// #endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
