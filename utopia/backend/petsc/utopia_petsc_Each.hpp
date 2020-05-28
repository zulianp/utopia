#ifndef UTOPIA_PETSC_EACH_HPP
#define UTOPIA_PETSC_EACH_HPP

//#include "utopia_Each.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Writable.hpp"
#include "utopia_petsc_RowView.hpp"
#include "utopia_petsc_Types.hpp"

#include <functional>

namespace utopia {

    template <int FILL_TYPE>
    class Each<PetscVector, 1, FILL_TYPE> {
    public:
        template <class Fun>
        inline static void apply_read(const PetscVector &v, Fun fun) {
            PetscErrorCode ierr;

            const auto r = range(v);
            const std::size_t r_begin = r.begin();
            const auto &impl = raw_type(v);

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

        template <class Fun>
        inline static void apply_write(PetscVector &v, Fun fun) {
            PetscErrorCode ierr;

            const auto r = range(v);
            const std::size_t r_begin = r.begin();
            const auto &impl = raw_type(v);

            PetscScalar *arr;

            ierr = VecGetArray(impl, &arr);
            assert(ierr == 0);

            For<>::apply(r_begin, r.end(), [&arr, &fun, r_begin](const std::size_t i) {
                auto idx = i - r_begin;
                arr[idx] = fun(i);
            });

            ierr = VecRestoreArray(impl, &arr);
            assert(ierr == 0);
            (void)ierr;
        }

        template <class Fun>
        inline static void apply_transform(PetscVector &in_out, Fun fun) {
            apply_transform(in_out, in_out, fun);
        }

        template <class Fun>
        inline static void apply_transform(const PetscVector &in, PetscVector &out, Fun fun) {
            PetscErrorCode ierr;

            const auto s = size(in);
            if (out.empty() || s != size(out)) {
                // out = local_zeros(local_size(in));
                out.zeros(layout(in));
            }

            assert(range(in) == range(out));

            const auto &impl_in = raw_type(in);
            auto &impl_out = raw_type(out);

            const auto r = range(out);
            const std::size_t r_begin = r.begin();

            if (in.is_alias(out)) {
                PetscScalar *arr;

                ierr = VecGetArray(impl_out, &arr);
                assert(ierr == 0);

                For<>::apply(r_begin, r.end(), [arr, &fun, r_begin](const std::size_t i) {
                    auto idx = i - r_begin;
                    arr[idx] = fun(i, arr[idx]);
                });

                ierr = VecRestoreArray(impl_out, &arr);
                assert(ierr == 0);

            } else {
                const PetscScalar *arr_in;
                PetscScalar *arr_out;

                ierr = VecGetArrayRead(impl_in, &arr_in);
                assert(ierr == 0);
                ierr = VecGetArray(impl_out, &arr_out);
                assert(ierr == 0);

                For<>::apply(r_begin, r.end(), [arr_in, arr_out, &fun, r_begin](const std::size_t i) {
                    auto idx = i - r_begin;
                    arr_out[idx] = fun(i, arr_in[idx]);
                });

                ierr = VecRestoreArrayRead(impl_in, &arr_in);
                assert(ierr == 0);
                ierr = VecRestoreArray(impl_out, &arr_out);
                assert(ierr == 0);
            }

            (void)ierr;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_EACH_HPP
