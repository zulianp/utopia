#ifndef UTOPIA_PETSC_PETSCCRSVIEW_HPP
#define UTOPIA_PETSC_PETSCCRSVIEW_HPP

#include <memory>

#include "utopia_ArrayView.hpp"

#include "petscmat.h"

namespace utopia {
    class PetscCrsView {
    public:
        PetscCrsView();
        PetscCrsView(Mat raw_mat);
        ~PetscCrsView();
        void set(Mat raw_mat);

        ArrayView<const PetscInt> row_ptr() const;
        ArrayView<const PetscInt> colidx() const;
        ArrayView<const PetscScalar> values() const;
        ArrayView<PetscScalar> values();

        PetscInt nnz() const;
        PetscInt cols() const;
        PetscInt rows() const;

        class Impl;

    private:
        std::shared_ptr<Impl> impl_;
    };

    template <>
    class Traits<PetscCrsView> {
    public:
        using Scalar = PetscScalar;
        using SizeType = PetscInt;
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_PETSCCRSVIEW_HPP
