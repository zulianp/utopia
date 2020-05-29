#ifndef UTOPIA_PETSC_TRILINOS_CONVERT_TENSOR_HPP
#define UTOPIA_PETSC_TRILINOS_CONVERT_TENSOR_HPP

#include "utopia_ConvertTensor.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"

namespace utopia {

    template <>
    class ConvertTensor<TpetraVector, PetscVector, 1, TRILINOS, PETSC> {
    public:
        static void apply(const TpetraVector &in, PetscVector &out);
    };

    template <>
    class ConvertTensor<PetscVector, TpetraVector, 1, PETSC, TRILINOS> {
    public:
        static void apply(const PetscVector &in, TpetraVector &out);
    };

    template <>
    class ConvertTensor<TpetraMatrix, PetscMatrix, 2, TRILINOS, PETSC> {
    public:
        static void apply(const TpetraMatrix &in, PetscMatrix &out);
    };

    template <>
    class ConvertTensor<PetscMatrix, TpetraMatrix, 2, PETSC, TRILINOS> {
    public:
        static void apply(const PetscMatrix &in, TpetraMatrix &out);
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_TRILINOS_CONVERT_TENSOR_HPP
