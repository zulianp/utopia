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
        void clear();
        bool empty() const;

        inline friend bool empty(const PetscCrsView &mat) { return mat.empty(); }

        ArrayView<const PetscInt> row_ptr() const;
        ArrayView<const PetscInt> colidx() const;
        ArrayView<const PetscScalar> values() const;
        ArrayView<PetscScalar> values();

        class RowView {
        public:
            inline RowView() : length(0), colidx_(nullptr), values_(nullptr) {}

            RowView(const PetscInt length, const PetscInt *colidx, PetscScalar *values)
                : length(length), colidx_(colidx), values_(values) {}

            const PetscInt length{0};

            inline PetscScalar value(const PetscInt idx) const {
                assert(idx < length);
                return values_[idx];
            }

            inline PetscScalar &value(const PetscInt idx) {
                assert(idx < length);
                return values_[idx];
            }

            inline PetscInt colidx(const PetscInt idx) const {
                assert(idx < length);
                return colidx_[idx];
            }

        private:
            const PetscInt *colidx_;
            PetscScalar *values_;
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

    void views_host(PetscMatrix &mat, PetscCrsView &d, PetscCrsView &o);
}  // namespace utopia

#endif  // UTOPIA_PETSC_PETSCCRSVIEW_HPP
