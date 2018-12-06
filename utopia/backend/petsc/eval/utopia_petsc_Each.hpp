#ifndef UTOPIA_PETSC_EACH_HPP
#define UTOPIA_PETSC_EACH_HPP

#include "utopia_Each.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {

	template<int FILL_TYPE>
	class Each<DVectord, 1, FILL_TYPE> {
	public:
		
		template<class Fun>
		inline static void apply_read(const DVectord &v, Fun fun)
		{
			PetscErrorCode ierr;

			const auto r = range(v);
			const std::size_t r_begin = r.begin();
			const auto &impl = raw_type(v);

			const PetscScalar *arr;

			ierr = VecGetArrayRead(impl, &arr); assert(ierr == 0);

			For<>::apply(
				r_begin,
				r.end(),
				[&arr, &fun, r_begin](const std::size_t i) {
					fun(r_begin + i, arr[i]);
				}
			);

			ierr = VecRestoreArrayRead(impl, &arr); assert(ierr == 0);
			(void) ierr;
		}

		template<class Fun>
		inline static void apply_write(DVectord &v, Fun fun)
		{
			PetscErrorCode ierr;

			const auto r = range(v);
			const std::size_t r_begin = r.begin();
			const auto &impl = raw_type(v);

			PetscScalar *arr;

			ierr = VecGetArray(impl, &arr); assert(ierr == 0);

			For<>::apply(
				r_begin,
				r.end(),
				[&arr, &fun, r_begin](const std::size_t i) {
					arr[i] = fun(r_begin + i);
				}
			);

			ierr = VecRestoreArray(impl, &arr); assert(ierr == 0);
			(void) ierr;
		}

		template<class Fun>
		inline static void apply_transform(const DVectord &in, DVectord &out, Fun fun)
		{
			PetscErrorCode ierr;

			const auto s = size(in);
			if(s != size(out)) {
				out = local_zeros(s);
			}

			const auto &impl_in = raw_type(in);
			auto &impl_out = raw_type(out);

			const auto r = range(out);
			const std::size_t r_begin = r.begin();

			if(impl_in == impl_out) {
				PetscScalar *arr;

				ierr = VecGetArray(impl_out, &arr); assert(ierr == 0);

				For<>::apply(
					r_begin,
					r.end(),
					[arr, &fun, r_begin](const std::size_t i) {
						arr[i] = fun(r_begin + i, arr[i]);
					}
				);

				ierr = VecRestoreArray(impl_out, &arr); assert(ierr == 0);

			} else {
				const PetscScalar *arr_in;
				PetscScalar *arr_out;

				ierr = VecGetArrayRead(impl_in, &arr_in); assert(ierr == 0);
				ierr = VecGetArray(impl_out, &arr_out);   assert(ierr == 0);

				For<>::apply(
					r_begin,
					r.end(),
					[arr_in, arr_out, &fun, r_begin](const std::size_t i) {
						arr_out[i] = fun(r_begin + i, arr_in[i]);
					}
				);

				ierr = VecRestoreArrayRead(impl_in, &arr_in); assert(ierr == 0);
				ierr = VecRestoreArray(impl_out, &arr_out);   assert(ierr == 0);
			}	

			(void) ierr;
		}
	};	
}

#endif //UTOPIA_PETSC_EACH_HPP
