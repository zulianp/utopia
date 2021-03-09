#ifndef UTOPIA_PETSC_ADDITIVE_CORRECTION_TRANSFER_HPP
#define UTOPIA_PETSC_ADDITIVE_CORRECTION_TRANSFER_HPP

#include "utopia_AdditiveCorrectionTransfer.hpp"

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Types.hpp"

// #include "utopia_petsc_BlockAdditiveCorrectionTransfer.hpp"

namespace utopia {

    template <>
    class AdditiveCorrectionTransfer<PetscMatrix, PetscVector, 1> final : public Transfer<PetscMatrix, PetscVector> {
    public:
        using Scalar = utopia::Traits<PetscVector>::Scalar;
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using IndexArray = utopia::Traits<PetscVector>::IndexArray;

        bool interpolate(const PetscVector &x_coarse, PetscVector &x_fine) const override;
        bool restrict(const PetscVector &x_fine, PetscVector &x_coarse) const override;
        bool restrict(const PetscMatrix &mat_fine, PetscMatrix &mat_coarse) const override;
        bool boolean_restrict_or(const PetscVector &, PetscVector &) override;
        bool project_down(const PetscVector &, PetscVector &) const override;

        bool project_down_positive_negative(const PetscVector &, const PetscVector &, PetscVector &) override;

        void init_memory() override;
        Scalar interpolation_inf_norm() const override;
        Scalar projection_inf_norm() const override;
        Scalar restriction_inf_norm() const override;

        void handle_equality_constraints(const PetscVector &) override;

        inline IndexArray &parent() { return parent_; }

        inline void set_size(const SizeType n_coarse_local, SizeType n_coarse_global) {
            n_coarse_local_ = n_coarse_local;
            n_coarse_global_ = n_coarse_global;
        }

    private:
        SizeType n_coarse_local_{-1}, n_coarse_global_{-1};
        IndexArray parent_;
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_ADDITIVE_CORRECTION_TRANSFER_HPP
