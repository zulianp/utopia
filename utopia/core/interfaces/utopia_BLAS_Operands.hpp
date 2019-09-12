#ifndef UTOPIA_BLAS_OPERANDS_HPP
#define UTOPIA_BLAS_OPERANDS_HPP

#include "utopia_Traits.hpp"

namespace utopia {

	template<class Tensor>
	class BLAS1Tensor {
	public:
		using Scalar   = typename utopia::Traits<Tensor>::Scalar;
		using SizeType = typename utopia::Traits<Tensor>::SizeType;
		
		virtual ~BLAS1Tensor() {}

		///<Scalar>SWAP - swap x and y
		virtual void swap(Tensor &x) = 0;

		///<Scalar>SCAL - x = a*x
		virtual void scale(const Scalar &a) = 0;

		///<Scalar>COPY - copy x into y (this)
		virtual void copy(const Tensor &x) = 0;

		///<Scalar>AXPY - y = a*x + y
		virtual void axpy(const Scalar &a, const Tensor &x) = 0;

		///<Scalar>DOT - dot product
		virtual Scalar dot(const Tensor &other) const = 0;

		///<Scalar>NRM2 - Euclidean norm
		virtual Scalar norm2() const = 0;

		///<Scalar>ASUM - sum of absolute values
		virtual Scalar asum() const = 0;

		///I<Scalar>AMAX - index of max abs value
		virtual SizeType amax() const = 0;

		//missing blas routines

		//<Scalar>ROTG - setup Givens rotation

		//<Scalar>ROTMG - setup modified Givens rotation

		//<Scalar>ROT - apply Givens rotation

		//<Scalar>ROTM - apply modified Givens rotation

		//<Scalar>SDOT - dot product with extended precision accumulation

		//<Scalar>ZNRM2 - Euclidean norm

	};

	template<class Matrix, class Vector>
	class BLAS2Matrix {
	public:
		using Scalar   = typename utopia::Traits<Matrix>::Scalar;
		using SizeType = typename utopia::Traits<Matrix>::SizeType;

		virtual ~BLAS2Matrix() {}

		/** @brief Original routine is
		  * <Scalar>GEMV - matrix vector multiply (y := alpha * A * x + beta * y)
		  * we separate split it in different methods since in general it is not 
		  * available in backends
		  */

		/// y = A * x
		virtual void multiply(const Vector &x, Vector &y) const
		{
			multiply(1.0, x, y);
		}

		/// y = alpha * A * x
		virtual void multiply(const Scalar &alpha, const Vector &x, Vector &y) const
		{
		    gemv(false, alpha, x, 0.0, y);
		}

		/// y = alpha * A * x
		virtual void transpose_multiply(const Vector &x, Vector &y) const
		{
		    gemv(false, 1.0, x, 0.0, y);
		}

		/// y = alpha * A^T * x
		virtual void transpose_multiply(const Scalar &alpha, const Vector &x, Vector &y) const
		{
		    gemv(true, alpha, x, 0.0, y);
		}

		/// y := alpha * A * x + beta * y
		virtual void multiply_add(const Scalar &alpha, const Vector &x, const Scalar &beta, Vector &y) const
		{
		    gemv(false, alpha, x, beta, y);
		}

		/// y := alpha * A' * x + beta * y
		virtual void transpose_multiply_add(const Scalar &alpha, const Vector &x, const Scalar &beta, Vector &y) const
		{
		    gemv(true, alpha, x, beta, y);
		}

		virtual void gemv(const bool transpose, const Scalar &alpha, const Vector &x, const Scalar &beta, Vector &y) const = 0;


		//missing blas routines

		// <Scalar>GBMV - banded matrix vector multiply

		// <Scalar>SYMV - symmetric matrix vector multiply

		// <Scalar>SBMV - symmetric banded matrix vector multiply

		// <Scalar>SPMV - symmetric packed matrix vector multiply

		// <Scalar>TRMV - triangular matrix vector multiply

		// <Scalar>TBMV - triangular banded matrix vector multiply

		// <Scalar>TPMV - triangular packed matrix vector multiply

		// <Scalar>TRSV - solving triangular matrix problems

		// <Scalar>TBSV - solving triangular banded matrix problems

		// <Scalar>TPSV - solving triangular packed matrix problems

		// <Scalar>GER - performs the rank 1 operation A := alpha*x*y' + A

		// <Scalar>SYR - performs the symmetric rank 1 operation A := alpha*x*x' + A

		// <Scalar>SPR - symmetric packed rank 1 operation A := alpha*x*x' + A

		// <Scalar>SYR2 - performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

		// <Scalar>SPR2 - performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

	};

	template<class Matrix>
	class BLAS3Matrix {
	public:
		using Scalar   = typename utopia::Traits<Matrix>::Scalar;
		using SizeType = typename utopia::Traits<Matrix>::SizeType;

		virtual ~BLAS3Matrix() {}

		virtual void transpose(Matrix &C) const = 0;

		virtual void multiply(const Matrix &B, Matrix &C) const
		{
			multiply(1.0, B, C);
		}

		/// C := alpha * A * B
		virtual void multiply(const Scalar &alpha, const Matrix &B, Matrix &C) const
		{
			multiply(false, alpha, false, B, C);
		}

		/// C := alpha * A * B
		virtual void transpose_multiply(const Matrix &B, Matrix &C) const
		{
			multiply(true, 1.0, false, B, C);
		}

		/// C := alpha * op(A) * op(B)
		virtual void multiply(
			const bool transpose_A,
			const bool transpose_B,
			const Matrix &B,
			Matrix &C) const
		{
			gemm(transpose_A, 1.0, transpose_B, B, 0.0, C);
		}

		/// C := alpha * op(A) * op(B)
		virtual void multiply(
			const bool transpose_A,
			const Scalar alpha,
			const bool transpose_B,
			const Matrix &B,
			Matrix &C) const
		{
			gemm(transpose_A, alpha, transpose_B, B, 0.0, C);
		}

		// <Scalar>GEMM - matrix matrix multiply  C := alpha*op( A )*op( B ) + beta*C
		virtual void gemm(
			const bool transpose_A,
			const Scalar alpha,
			const bool transpose_B,
			const Matrix &B,
			const Scalar beta,
			Matrix &C)  const = 0;

		//missing blas routines

		// <Scalar>SYMM - symmetric matrix matrix multiply

		// <Scalar>SYRK - symmetric rank-k update to a matrix

		// <Scalar>SYR2K - symmetric rank-2k update to a matrix

		// <Scalar>TRMM - triangular matrix matrix multiply

		// <Scalar>TRSM - solving triangular matrix with multiple right hand sides
	};
}

#endif //UTOPIA_BLAS_OPERANDS_HPP
