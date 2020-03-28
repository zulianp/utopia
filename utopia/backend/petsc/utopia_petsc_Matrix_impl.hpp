#ifndef UTOPIA_PETSC_MATRIX_IMPL_HPP
#define UTOPIA_PETSC_MATRIX_IMPL_HPP

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Each.hpp"

//Transform calue mpiaij
//1)  MatMPIAIJGetSeqAIJ (not MatMPIAIJGetLocalMat)
// PETSC_EXTERN PetscErrorCode MatMPIAIJGetSeqAIJ(Mat,Mat*,Mat*,const PetscInt*[]);

#define UNROLL_FACTOR 5


namespace utopia {

	namespace internal {

		template<class F>
		void read_petsc_seqaij_impl(Mat &mat, F op)
		{
			PetscInt n;
			const PetscInt *ia;
			const PetscInt *ja;
			PetscBool done;
			PetscErrorCode err = 0;

			err = MatGetRowIJ(mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
			assert(done == PETSC_TRUE);

			if(!done) {
			    std::cerr << "PetscMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what was asked." << std::endl;
			    abort();
			}

			PetscScalar *array;
			MatSeqAIJGetArray(mat, &array);

			for(PetscInt i = 0; i < n; ++i) {
				const PetscInt row_end = ia[i+1];

				#pragma clang loop unroll_count(UNROLL_FACTOR)
				#pragma GCC unroll 5
				for(PetscInt k = ia[i]; k < row_end; ++k) {
			    	op(i, ja[k], array[k]);
			    }
			}

			MatSeqAIJRestoreArray(mat, &array);
			err = MatRestoreRowIJ(mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
		}

		template<class F>
		void read_reverse_petsc_seqaij_impl(Mat &mat, F op)
		{
			PetscInt n;
			const PetscInt *ia;
			const PetscInt *ja;
			PetscBool done;
			PetscErrorCode err = 0;

			err = MatGetRowIJ(mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
			assert(done == PETSC_TRUE);

			if(!done) {
			    std::cerr << "PetscMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what was asked." << std::endl;
			    abort();
			}

			PetscScalar *array;
			MatSeqAIJGetArray(mat, &array);

			for(PetscInt i = n-1; i >= 0; --i) {
				const PetscInt row_end = ia[i+1];

				#pragma clang loop unroll_count(UNROLL_FACTOR)
				#pragma GCC unroll 5
				for(PetscInt k = ia[i]; k < row_end; ++k) {
			    	op(i, ja[k], array[k]);
			    }
			}

			MatSeqAIJRestoreArray(mat, &array);
			err = MatRestoreRowIJ(mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
		}

		template<class F>
		void read_petsc_mpiaij_impl(Mat &mat, F op)
		{
			// PetscErrorCode err = 0;

			// const PetscInt* cols;
			// Mat d, o;
			// err =  MatMPIAIJGetSeqAIJ(mat, &d, &o, &cols); assert(err == 0);

			assert(false && "IMPLEMENT ME");
		}

		template<class F>
		void transform_petsc_impl(Mat &mat, F op)
		{
			PetscInt n;
			const PetscInt *ia;
			// const PetscInt *ja;
			PetscBool done;
			PetscErrorCode err = 0;

			err = MatGetRowIJ(mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, nullptr, &done); assert(err == 0);
			assert(done == PETSC_TRUE);

			if(!done) {
			    std::cerr << "PetscMatrix::transform_petsc_impl(const Op &op): MatGetRowIJ failed to provide what was asked." << std::endl;
			    abort();
			}

			PetscScalar *array;
			MatSeqAIJGetArray(mat, &array);

			//FIXME is there a better way to get the total number of values???
			const PetscInt n_values = ia[n] - ia[0];

			for(PetscInt i = 0; i < n_values; ++i) {
			    array[i] = op(array[i]);
			}

			MatSeqAIJRestoreArray(mat, &array);
			err = MatRestoreRowIJ(mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, nullptr, &done); assert(err == 0);
		}
	}

	template<class F>
	void PetscMatrix::transform_values_mpiaij(F op)
	{
		PetscErrorCode err = 0;

		const PetscInt* cols;
		Mat d, o;
		err =  MatMPIAIJGetSeqAIJ(raw_type(), &d, &o, &cols); assert(err == 0);

		internal::transform_petsc_impl(d, op);
		internal::transform_petsc_impl(o, op);

		//DOUBT no restore?
		// err =  MatMPIAIJRestoreSeqAIJ(raw_type(), &d, &o, &cols); assert(err == 0);
	}

	template<class F>
	void PetscMatrix::transform_values_seqaij(F op)
	{
		internal::transform_petsc_impl(raw_type(), op);
	}

	template<class F>
	void PetscMatrix::transform_values(F op)
	{
	    if(has_type(MATSEQAIJ)) {
	        transform_values_seqaij(op);
	    } else if(has_type(MATMPIAIJ)) {
	    	transform_values_mpiaij(op);
	    } else {
	        each_transform(*this, [op](const SizeType &, const SizeType &, const Scalar value) -> Scalar {
	            return op(value);
	        });
	    }
	}

	template<class F>
	void PetscMatrix::transform_ijv_seqaij(F op)
	{
		PetscInt n;
		const PetscInt *ia;
		const PetscInt *ja;
		PetscBool done;
		PetscErrorCode err = 0;

		err = MatGetRowIJ(raw_type(), 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
		assert(done == PETSC_TRUE);

		if(!done) {
		    std::cerr << "PetscMatrix::transform_values_seqaij(const Op &op): MatGetRowIJ failed to provide what was asked." << std::endl;
		    abort();
		}

		PetscScalar *array;
		MatSeqAIJGetArray(raw_type(), &array);

		auto ra = this->row_range();
		PetscInt n_local_rows = ra.extent();
		for(PetscInt r = 0; r < n_local_rows; ++r) {
			const PetscInt row_begin = ia[r];
			const PetscInt row_end   = ia[r+1];

			const PetscInt n_values = row_end - row_begin;
			const PetscInt row_global = ra.begin() + r;

			for(PetscInt i = row_begin; i < row_end; ++i) {
			    array[i] = op(row_global, ja[i], array[i]);
			}
		}

		MatSeqAIJRestoreArray(raw_type(), &array);
		err = MatRestoreRowIJ(raw_type(), 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
	}

	//wait for petsc version
	// template<class F>
	// void PetscMatrix::read_ijv_seqaij(F op)
	// {
	// 	PetscInt n;
	// 	const PetscInt *ia;
	// 	const PetscInt *ja;
	// 	PetscBool done;
	// 	PetscErrorCode err = 0;

	// 	err = MatGetRowIJ(raw_type(), 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
	// 	assert(done == PETSC_TRUE);

	// 	if(!done) {
	// 	    std::cerr << "PetscMatrix::transform_values_seqaij(const Op &op): MatGetRowIJ failed to provide what was asked." << std::endl;
	// 	    abort();
	// 	}

	// 	const PetscScalar *array;
	// 	MatSeqAIJGetArrayRead(raw_type(), &array);

	// 	auto ra = this->row_range();
	// 	PetscInt n_local_rows = ra.extent();
	// 	for(PetscInt r = 0; r < n_local_rows; ++r) {
	// 		const PetscInt row_begin = ia[r];
	// 		const PetscInt row_end   = ia[r+1];

	// 		const PetscInt n_values = row_end - row_begin;
	// 		const PetscInt row_global = ra.begin() + r;

	// 		for(PetscInt i = row_begin; i < row_end; ++i) {
	// 		    op(row_global, ja[i], array[i]);
	// 		}
	// 	}

	// 	MatSeqAIJRestoreArrayRead(raw_type(), &array);
	// 	err = MatRestoreRowIJ(raw_type(), 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done); assert(err == 0);
	// }

	template<class F>
	void PetscMatrix::transform_ijv(F op)
	{
	    if(has_type(MATSEQAIJ)) {
	        transform_ijv_seqaij(op);
	    } else {
	        each_transform(*this, [op](const SizeType &i, const SizeType &j, const Scalar value) -> Scalar {
	            return op(i, j, value);
	        });
	    }
	}

	template<class Op>
	void PetscMatrix::op_transform(const Op &op)
	{
	    transform_values([op](const Scalar &value) -> Scalar {
	        return op.template apply(value);
	    });
	}

	template<class Op>
	void PetscMatrix::read(Op op)
	{
		if(has_type(MATSEQAIJ)) {
			internal::read_petsc_seqaij_impl(raw_type(), op);
		} else {
			assert(false);
		}
	}

	template<class Op>
	void PetscMatrix::read_reverse(Op op)
	{
		if(has_type(MATSEQAIJ)) {
			internal::read_reverse_petsc_seqaij_impl(raw_type(), op);
		} else {
			assert(false);
		}
	}

	template<class Operation>
	inline static PetscMatrix::Scalar generic_local_reduce(const PetscMatrix &m, const PetscMatrix::Scalar &init_value, const Operation &op)
	{
	    using Scalar = PetscMatrix::Scalar;

	    Scalar x = init_value;
	    const Scalar * values;
	    const PetscInt * cols;

	    PetscInt r_begin, r_end;
	    PetscInt n_values = 0;

	    PetscInt local_r, local_c;

	    MatGetLocalSize(m.raw_type(), &local_r, &local_c);
	    MatGetOwnershipRange(m.raw_type(), &r_begin, &r_end);

	    for(PetscInt row = r_begin; row < r_end; ++row) {

	        MatGetRow(m.raw_type(), row, &n_values, &cols, &values);

	        if(n_values < local_c) {
	            x = op.template apply<Scalar>(x, 0.);
	        }

	        for(PetscInt i = 0; i < n_values; ++i) {
	            x = op.template apply<Scalar>(x, values[i]);
	        }

	        MatRestoreRow(m.raw_type(), row, &n_values, &cols, &values);
	    }

	    return x;
	}

	template<class Operation>
	inline static void reduce_rows(PetscVector &result,
	 const PetscMatrix &mat,
	 const PetscMatrix::Scalar &init_value,
	 const Operation &op
	 )
	{
	    using Scalar = PetscMatrix::Scalar;

	    assert(!result.is_null());

	    const Scalar * values;
	    const PetscInt * cols;

	    PetscInt r_begin, r_end;
	    PetscInt n_values = 0;

	    PetscInt global_r, global_c, local_r, local_c;
	    MatGetSize(mat.raw_type(), &global_r, &global_c);
	    MatGetLocalSize(mat.raw_type(), &local_r, &local_c);

	    MatGetOwnershipRange(mat.raw_type(), &r_begin, &r_end);

	    result.write_lock(utopia::LOCAL);

	    for(PetscInt row = r_begin; row < r_end; ++row) {
	        MatGetRow(mat.raw_type(), row, &n_values, &cols, &values);

	        Scalar x = init_value;
	        for(PetscInt i = 0; i < n_values; ++i) {
	            x = op.template apply<Scalar>(x, values[i]);
	        }

	        if(n_values < global_c) {
	            x = op.template apply<Scalar>(x, 0.);
	        }

	        MatRestoreRow(mat.raw_type(), row, &n_values, &cols, &values);
	        VecSetValues(result.raw_type(), 1, &row, &x, INSERT_VALUES);
	    }

	    result.write_unlock(utopia::LOCAL);
	}

}

#undef UNROLL_FACTOR //clean-up
#endif //UTOPIA_PETSC_MATRIX_IMPL_HPP

