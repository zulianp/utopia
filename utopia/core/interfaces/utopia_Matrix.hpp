#ifndef UTOPIA_MATRIX_HPP
#define UTOPIA_MATRIX_HPP

#include "utopia_Enums.hpp"

namespace utopia {

	template<typename Scalar_, typename SizeType_>
	class MatrixBase {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~MatrixBase() {}

		//locks
		virtual void read_lock() 			 	= 0;
		virtual void write_lock(WriteMode mode) = 0;

		virtual void read_unlock()  			  = 0;
		virtual void write_unlock(WriteMode mode) = 0;

		//basic mutators
		virtual void set(const SizeType &i, const SizeType &j, const Scalar &value) = 0;
		virtual void add(const SizeType &i, const SizeType &j, const Scalar &value) = 0;

		//print function
		virtual void describe() const = 0;

		//utility functions
		virtual bool empty() const = 0;
		virtual bool clear() const = 0;
	};

	template<typename Scalar_, typename SizeType_>
	class Matrix : public MatrixBase<Scalar_, SizeType_> {
	public:
		virtual ~Matrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class DenseMatrix : public Matrix<Scalar_, SizeType_> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual Scalar get(const SizeType &i, const SizeType &j) = 0;

		virtual ~DenseMatrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class SparseMatrix : public Matrix<Scalar_, SizeType_> {
	public:
		virtual ~SparseMatrix() {}
	};

	//parallel types, collective operations
	template<typename Scalar_, typename SizeType_>
	class DistributedMatrix : public MatrixBase<Scalar_, SizeType_> {
	public:
		using Scalar = Scalar_;
		using SizeType = SizeType_;

		//basic collective mutators allowing to write on other processes (e.g. for FE assembly)
		virtual void c_set(const SizeType &i, const SizeType &j, const Scalar &value) = 0;
		virtual void c_add(const SizeType &i, const SizeType &j, const Scalar &value) = 0;

		virtual ~DistributedMatrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class DistributedDenseMatrix : public DistributedMatrix<Scalar_, SizeType_> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual Scalar get(const SizeType &i, const SizeType &j) = 0;

		virtual ~DistributedDenseMatrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class DistributedSparseMatrix : public DistributedMatrix<Scalar_, SizeType_> {
	public:
		virtual ~DistributedSparseMatrix() {}
	};

}

#endif //UTOPIA_MATRIX_HPP
