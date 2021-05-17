#ifndef UTOPIA_PETSC_PETSCCRSVIEW_HPP
#define UTOPIA_PETSC_PETSCCRSVIEW_HPP

#include <memory>

#include "utopia_ArrayView.hpp"
#include "utopia_Describable.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

#include "petscmat.h"

namespace utopia {
    class PetscCrsView : public Describable {
    public:
        PetscCrsView();
        PetscCrsView(Mat raw_mat);
        ~PetscCrsView();
        void set(Mat raw_mat);

        ArrayView<const PetscInt> row_ptr() const;
        ArrayView<const PetscInt> colidx() const;
        ArrayView<const PetscScalar> values() const;
        ArrayView<PetscScalar> values();

        class RowView {
        public:
            RowView(const PetscInt length, const PetscInt *colidx, const PetscScalar *values)
                : length(length), colidx_(colidx), values_(values) {}

            const PetscInt length{0};

            inline PetscScalar value(const PetscInt idx) const {
                assert(idx < length);
                return values_[idx];
            }

            inline PetscInt colidx(const PetscInt idx) const {
                assert(idx < length);
                return colidx_[idx];
            }

        private:
            const PetscInt *colidx_;
            const PetscScalar *values_;
        };

        RowView row(const PetscInt row) const;

        PetscInt nnz() const;
        PetscInt cols() const;
        PetscInt rows() const;

        class Impl;

        void describe(std::ostream &os = std::cout) const override;

    private:
        std::shared_ptr<Impl> impl_;
    };

    template <>
    class Traits<PetscCrsView> {
    public:
        using Scalar = PetscScalar;
        using SizeType = PetscInt;
    };

    PetscCrsView crs_view(PetscMatrix &mat);
    PetscCrsView crs_view(const PetscMatrix &mat);
}  // namespace utopia

#endif  // UTOPIA_PETSC_PETSCCRSVIEW_HPP
