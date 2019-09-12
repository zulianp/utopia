#ifndef UTOPIA_MATRIX_HPP
#define UTOPIA_MATRIX_HPP

#include "utopia_Enums.hpp"
#include "utopia_Range.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_DistributedObject.hpp"
#include "utopia_Size.hpp"

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
		virtual void clear() = 0;


		virtual SizeType rows() const = 0;
		virtual SizeType cols() const = 0;
	};

	template<typename Scalar_, typename SizeType_>
	class Matrix : public MatrixBase<Scalar_, SizeType_> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual Size size() const = 0;
		virtual ~Matrix() {}

		//face methods for treaing matrix and distrubed-matrix in the same way
		void c_set(const SizeType &i, const SizeType &j, const Scalar &value)
		{
			this->set(i, j, value);
		}

		void c_add(const SizeType &i, const SizeType &j, const Scalar &value)
		{
			this->add(i, j, value);
		}

		inline void row_range(Range &r) const { r.set(0, this->rows()); }
		inline void col_range(Range &r) const { r.set(0, this->cols()); }

		inline void row_layout(Layout<SizeType, 1> &l) const
		{
			l.local_size(0) = local_rows();
			l.global_size(0) = this->rows();
		}

		inline void col_layout(Layout<SizeType, 1> &l) const
		{
			l.local_size(0) = local_cols();
			l.global_size(0) = this->cols();
		}

		inline void layout(Layout<SizeType, 2> &l) const
		{
			l.local_size(0) = local_rows();
			l.local_size(1) = local_cols();

			l.global_size(0) = this->rows();
			l.global_size(1) = this->cols();
		}

		inline SizeType local_rows() const { return this->rows(); }
		inline SizeType local_cols() const { return this->cols(); }

		inline Size local_size() const { return {local_rows(), local_cols() }; }

	};

	template<typename Scalar_, typename SizeType_>
	class DenseMatrix : public Matrix<Scalar_, SizeType_> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;
		using Super    = utopia::Matrix<Scalar_, SizeType_>;
		using Super::set;

		virtual Scalar get(const SizeType &i, const SizeType &j) const = 0;
		virtual void set(const Scalar &val) = 0;

		virtual ~DenseMatrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class SparseMatrix : public Matrix<Scalar_, SizeType_> {
	public:
		virtual ~SparseMatrix() {}
	};

	//parallel types, collective operations
	template<typename Scalar_, typename SizeType_>
	class DistributedMatrix : public MatrixBase<Scalar_, SizeType_>, public DistributedObject {
	public:
		using Scalar = Scalar_;
		using SizeType = SizeType_;

		//basic collective mutators allowing to write on other processes (e.g. for FE assembly)
		virtual void c_set(const SizeType &i, const SizeType &j, const Scalar &value) = 0;
		virtual void c_add(const SizeType &i, const SizeType &j, const Scalar &value) = 0;

		virtual void row_range(Range &r) const = 0;
		virtual void col_range(Range &r) const = 0;

		virtual void row_layout(Layout<SizeType, 1> &l) const
		{
			l.local_size(0) = local_rows();
			l.global_size(0) = this->rows();
		}

		virtual void col_layout(Layout<SizeType, 1> &l) const
		{
			l.local_size(0) = local_cols();
			l.global_size(0) = this->cols();
		}

		virtual void layout(Layout<SizeType, 2> &l) const
		{
			l.local_size(0) = local_rows();
			l.local_size(1) = local_cols();

			l.global_size(0) = this->rows();
			l.global_size(1) = this->cols();
		}

		virtual SizeType local_rows() const = 0;
		virtual SizeType local_cols() const = 0;
		virtual Size local_size() const { return {local_rows(), local_cols() }; }

		virtual ~DistributedMatrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class DistributedDenseMatrix : public DistributedMatrix<Scalar_, SizeType_> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual Scalar get(const SizeType &i, const SizeType &j) const = 0;

		virtual ~DistributedDenseMatrix() {}
	};

	template<typename Scalar_, typename SizeType_>
	class DistributedSparseMatrix : public DistributedMatrix<Scalar_, SizeType_> {
	public:
		virtual ~DistributedSparseMatrix() {}
	};

}

#endif //UTOPIA_MATRIX_HPP
