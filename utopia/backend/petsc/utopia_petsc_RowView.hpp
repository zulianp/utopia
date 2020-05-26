#ifndef UTOPIA_PETSC_ROW_VIEW_HPP
#define UTOPIA_PETSC_ROW_VIEW_HPP

#include "utopia_RowView.hpp"

namespace utopia {

    template <class Tensor, int FILL_TYPE>
    class RowView<Tensor, 2, FILL_TYPE, utopia::PETSC> {
    public:
        // FIXME move petsc specific stuff to backend

        inline RowView(Tensor &t, const PetscInt row) : t_(t), row_(row) {
            MatGetRow(raw_type(t_), row_, &n_values_, &cols_, &values_);
        }

        inline ~RowView() { MatRestoreRow(raw_type(t_), row_, &n_values_, &cols_, &values_); }

        inline PetscInt n_values() const { return n_values_; }

        inline PetscInt col(const PetscInt index) const {
            assert(index < n_values_);
            return cols_[index];
        }

        inline PetscScalar get(const PetscInt index) const {
            assert(index < n_values_);
            return values_[index];
        }

    private:
        Tensor &t_;
        PetscInt row_;
        PetscInt n_values_{0};
        const PetscInt *cols_{nullptr};
        const PetscScalar *values_{nullptr};
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_ROW_VIEW_HPP
