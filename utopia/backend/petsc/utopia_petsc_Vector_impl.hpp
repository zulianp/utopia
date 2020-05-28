#ifndef UTOPIA_PETSC_VECTOR_IMPL_HPP
#define UTOPIA_PETSC_VECTOR_IMPL_HPP

#include "utopia_petsc_Vector.hpp"

namespace utopia {

    template <class Operation>
    void PetscVector::op_transform(const Operation &op) {
        transform_values([op](const Scalar &val) -> Scalar { return op.template apply<Scalar>(val); });
    }

    template <class F>
    void PetscVector::transform_values(F op) {
        read_and_write_lock(LOCAL);

        auto values = writeable_->data;
        SizeType n = writeable_->range_end - writeable_->range_begin;

        for (SizeType i = 0; i < n; ++i) {
            auto &v = values[i];
            v = op(v);
        }

        read_and_write_unlock(LOCAL);
    }

    template <typename Scalar, class Operation>
    static void element_wise_generic(const Scalar &x, const Operation &op, PetscVector &result) {
        result.transform_values([op, x](const Scalar &val) -> Scalar { return op.template apply<Scalar>(x, val); });
    }

    template <class Fun>
    void PetscVector::read(Fun fun) const {
        PetscErrorCode ierr;

        const auto r = this->range();
        const std::size_t r_begin = r.begin();
        const auto &impl = this->raw_type();

        const PetscScalar *arr;

        ierr = VecGetArrayRead(impl, &arr);
        assert(ierr == 0);

        For<>::apply(r_begin, r.end(), [&arr, &fun, r_begin](const std::size_t i) {
            auto idx = i - r_begin;
            fun(i, arr[idx]);
        });

        ierr = VecRestoreArrayRead(impl, &arr);
        assert(ierr == 0);
        (void)ierr;
    }

}  // namespace utopia

#endif  // UTOPIA_PETSC_VECTOR_IMPL_HPP
