#ifndef utopia_utopia_BLASBACKEND_HPP
#define utopia_utopia_BLASBACKEND_HPP

#include "utopia_blas_Matrix.hpp"
#include "utopia_blas_Traits.hpp"

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Core.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_ScalarBackend.hpp"

#include <vector>
#include <memory>
#include <iostream>
#include <numeric>
#include <algorithm>

namespace utopia {
	
	class BLASBackend : public ScalarBackend<double> {
	public:
		typedef double Scalar;
		typedef utopia::Matrix<Scalar> Matrix;
		typedef std::vector<Scalar> Vector;
		typedef std::vector<Scalar>::size_type SizeType;
		
		using ScalarBackend<Scalar>::apply_binary;
		using ScalarBackend<Scalar>::axpy;
		
		void static daxpy_wrapper(const int *n, const Scalar *alpha, const Scalar *x, const int *incx, Scalar *y, const int *incy);
		void static dscal_wrapper(const int *n, const Scalar *sa, const Scalar *sx, const int *incx);
		
		template<class Tensor>
		inline static void read_lock(const Tensor &) {}
		
		template<class Tensor>
		inline static void read_unlock(const Tensor &) {}
		
		template<class Tensor>
		inline static void write_lock(const Tensor &) {}
		
		template<class Tensor>
		inline static void write_unlock(const Tensor &) {}
		
		inline static void write_lock(CRSMatrix<Scalar> &mat) {
			mat.assembly_begin();
		}
		
		inline static void write_unlock(CRSMatrix<Scalar> &mat) {
			mat.assembly_end();
		}
		
		inline static void write_lock(CCSMatrix<Scalar> &mat) {
			mat.assembly_begin();
		}
		
		inline static void write_unlock(CCSMatrix<Scalar> &mat) {
			mat.assembly_end();
		}

		template<class Tensor>
		inline static void read_and_write_lock(const Tensor &) {}
		
		template<class Tensor>
		inline static void read_and_write_unlock(const Tensor &) {}
		
		inline static Range range(const Vector &v) {
			return Range(0, v.size());
		}
		
		template<class MatrixT>
		inline static Range row_range(const MatrixT &m) {
			return Range(0, m.rows());
		}
		
		template<class MatrixT>
		inline static Range col_range(const MatrixT &m) {
			return Range(0, m.cols());
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//[queries]
		template<class Tensor>
		static void local_size(const Tensor &t, Size &result)
		{
			return size(t, result);
		}
		
		static void size(const Vector &v, Size &size);
		static void size(const Matrix &m, Size &size);
		static void size(const CRSMatrix<Scalar> &m, Size &size);
		static void size(const CCSMatrix<Scalar> &m, Size &size);
		
		//[builders]
		
		template<class Tensor>
		static void build(Tensor &t, const Size &s, const Resize &)
		{
			resize(t, s);
		}
		
		static void build(Matrix &m, const Size &size, const Identity &);
		static void build(CRSMatrix<Scalar> &m, const Size &size, const Identity &);
		static void build(CRSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz);
		static void build(CCSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz);
		static void build(CRSMatrix<Scalar> &m, const Size &size, const Zeros & /*values*/);
		static void build(Matrix &m, const Size &size, const Values<Scalar> &values);
		static void build(Matrix &m, const Size &size, const Zeros & /*values*/);
		static void build(Vector &m, const Size &size, const Zeros & /*values*/);
		static void build(Vector &m, const Size &size, const LocalZeros & /*values*/);
		static void build(Vector &v, const Size &size, const Values<Scalar> &values);
		
		//[accessors]
		inline static Scalar get(const Vector &vec, const SizeType index)
		{
			return vec[index];
		}
		
		//[mutators]
		template<class Left>
		static void apply(Left &t, const Resize &, const Size &s)
		{
			resize(t, s);
		}
		
		static void resize(Matrix &mat, const Size &size);
		static void resize(Vector &vec, const Size &size);
		
		inline static void set(Vector &vec, const SizeType index, const Scalar value)
		{
			assert(index < vec.size());
			vec[index] = value;
		}
		
		template<class MatrixT>
		inline static void set(
							   MatrixT &mat,
							   const SizeType row,
							   const SizeType col,
							   const Scalar value)
		{
			assert(row < mat.rows() && col < mat.cols());
			mat.set(row, col, value);
		}
		
		template<typename Ordinal>
		inline static void set(
							   Vector &v,
							   const std::vector<Ordinal> &indices,
							   const std::vector<Scalar> &values)
		{
			assert(indices.size() == values.size());
			
			for (typename std::vector<Ordinal>::size_type i = 0; i < indices.size(); ++i)
				set(v, indices[i], values[i]);
		}


		static void set(Vector &v, Scalar value)
		{
			std::fill(v.begin(), v.end(), value);
		}
		
		// template<typename Ordinal>
		// inline static void set(
		// 					   Matrix &m,
		// 					   const std::vector<Ordinal> &rows,
		// 					   const std::vector<Ordinal> &columns,
		// 					   const std::vector<Scalar> &values)
		// {
		// 	assert(rows.size() == columns.size() && rows.size() == values.size());
			
		// 	for (typename std::vector<Ordinal>::size_type i = 0; i < rows.size(); ++i)
		// 		set(m, rows[i], columns[i], values[i]);
		// }

		template<typename Ordinal>
		void set_matrix(Matrix &m,
							   const std::vector<Ordinal> &rows,
							   const std::vector<Ordinal> &cols,
							   const std::vector<Scalar> &values)
		{

			const auto n_rows = rows.size();
			const auto n_cols = cols.size();

			for(std::size_t i = 0; i < n_rows; ++i) {
				for(std::size_t j = 0; j < n_cols; ++j) {
					set(m, rows[i], cols[j], values[i * n_rows + j]);
				}
			}
		}
		
		template<typename Ordinal>
		inline static void set(
							   CRSMatrix<Scalar> &m,
							   const std::vector<Ordinal> &rows,
							   const std::vector<Ordinal> &columns,
							   const std::vector<Scalar> &values)
		{
			for (typename std::vector<Ordinal>::size_type i = 0; i < rows.size(); ++i) {
				set(m, rows[i], columns[i], values[i]);
			}
		}
		
		template<class MatrixT>
		inline static Scalar get(const MatrixT &mat, const SizeType row, const SizeType col) {
			assert(row < mat.rows() && col < mat.cols());
			return mat.get(row, col);
		}
		
		inline static void add(Vector &vec, const SizeType index, const Scalar value)
		{
			assert(index < vec.size());
			vec[index] += value;
		}
		
		template<class MatrixT>
		inline static void add(MatrixT &mat, const SizeType row, const SizeType col, const Scalar value) {
			assert(row < mat.rows() && col < mat.cols());
			mat.set(row, col, mat.get(row, col) + value);
		}
		
		//[boolean operations]
		template<typename T, class Comparator>
		inline static bool compare(const std::vector<T> &left, const std::vector<T> &right, const Comparator &comp) {
			assert(left.size() == right.size());
			
			if (left.size() != right.size()) return false;
			
			bool result = true;
			for (SizeType i = 0; i < left.size(); ++i) {
				result &= comp(left[i], right[i]);
			}
			
			return result;
		}
		
		template<class Comparator>
		inline static bool compare(const Matrix &left, const Matrix &right, const Comparator &comp) {
			assert(left.rows() == right.rows());
			assert(left.cols() == right.cols());
			
			if (left.rows() != right.rows() || left.cols() != right.cols()) return false;
			
			return compare(left.entries(), right.entries(), comp);
		}
		
		template<class Comparator>
		inline static bool compare(const CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right, const Comparator &comp) {
			assert(left.rows() == right.rows());
			assert(left.cols() == right.cols());
			
			if (left.rows() != right.rows() || left.cols() != right.cols()) return false;
			
			for (SizeType r = 0; r != left.rows(); ++r) {
				const SizeType boundaryl = left.rowptr()[r + 1];
				const SizeType boundaryr = right.rowptr()[r + 1];
				
				for (SizeType kl = left.rowptr()[r], kr = right.rowptr()[r];
					 kl != boundaryl && kr != boundaryr;
					 /* increment inside */
					 )
				{
					const SizeType cl = left.colindex()[kl];
					const SizeType cr = right.colindex()[kr];
					const Scalar vl = left.at(kl);
					const Scalar vr = right.at(kr);
					
					if(cl == cr) {
						if (!comp(vl, vr)) {
							return false;
						} else {
							++kl;
							++kr;
						}
					}
					else if(cl < cr)
					{
						if(!comp(vl, Scalar(0))) {
							return false;
						}
						++kl;
					}
					else if(cl > cr) {
						if(!comp(vr, Scalar(0))) {
							return false;
						}
						++kr;
					}
					
					if(kl == boundaryl && kr != boundaryr) {
						for(; kl != boundaryl; ++kl) {
							if(!comp(left.at(kl), Scalar(0))) {
								return false;
							}
						}
						
					} else if(kl != boundaryl && kr == boundaryr) {
						for(; kr != boundaryr; ++kr) {
							if(!comp(right.at(kr), Scalar(0))) {
								return false;
							}
						}
					}
				}
			}
			
			return true;
		}
		
		//[utilities]
		inline static bool handle_error(int err) {
			ASSERT(err == 0);
			return err == 0;
		}
		
		inline static Scalar *ptr(Matrix &mat)
		{
			return mat.ptr();
		}
		
		inline static Scalar *ptr(Vector &v)
		{
			return &v[0];
		}
		
		inline static const Scalar *ptr(const Matrix &mat)
		{
			return mat.ptr();
		}
		
		inline static const Scalar *ptr(const Vector &v)
		{
			return &v[0];
		}
		
		inline static const Scalar *ptr(const CRSMatrix<Scalar> &mat)
		{
			return &mat.entries()[0];
		}
		
		template<typename T>
		inline static const T *ptr(const std::vector<T> &v) {
			return &v[0];
		}
		
		//[assignments]
		//vector
		static void assign(Vector &left, Vector &&right);
		static void assign(Vector &left, const Vector &right);
		
		
		//matrix
		static void assign(Matrix &left, Matrix &&right) ;
		static void assign(Matrix &left, const CRSMatrix<Scalar> &right);
		static void assign(Matrix &left, const Matrix &right);
		static void assign(CCSMatrix<Scalar> &left, const CCSMatrix<Scalar> &right);
		static void assign(CCSMatrix<Scalar> &left, CCSMatrix<Scalar> &&right);
		static void assign(CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right);
		static void assign(CRSMatrix<Scalar> &left, CRSMatrix<Scalar> &&right);
		static void assign(CCSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right);

        static void assign_transposed(Matrix &left, const CRSMatrix<Scalar> &right);
        static void assign_transposed(Matrix &left, const Matrix &right);
        static void assign_transposed(CCSMatrix<Scalar> &left, const CCSMatrix<Scalar> &right);
        static void assign_transposed(CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right);
        static void assign_transposed(CCSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right);
		
		//[blas 1]
		static Scalar trace(const Matrix &in);
		static Scalar trace(const CCSMatrix<Scalar>  &in);
		static Scalar dot(const Vector &left, const Vector &right);
		static Scalar dot(const Matrix &left, const Matrix &right);
		static Scalar norm2(const Vector &vector);
		static Scalar norm2(const Matrix &mat);
		static Scalar norm2(const CRSMatrix<Scalar> &mat);
		static Scalar norm2(const CCSMatrix<Scalar> &mat);
		static Scalar norm_infty(const Vector &vector);
		static Scalar reduce(const Vector &vec, const Plus &);
		static Scalar reduce(const Matrix &m,   const Plus &);
		
		//axpy
		static void axpy(Vector &y, const Scalar &alpha, const Vector &x);
		static void axpy(Matrix &y, const Scalar &alpha, const Matrix &x);


		//scale
		static void scale(Vector &result, const Scalar scale_factor);
		static void scale(Matrix &result, const Scalar scale_factor);
        //TODO
		// static void axpy(CCSMatrix<Scalar> &y, const Scalar &alpha, const CCSMatrix<Scalar> &x);
		// static void axpy(CRSMatrix<Scalar> &y, const Scalar &alpha, const CRSMatrix<Scalar> &x);
		
		//[blas 2]

		inline static void multiply(Vector &y,
						 const bool transpose_A,
						 const Matrix &A,
						 const bool transpose_X,
						 const Vector &x)
		{
			assert(!transpose_X);
			(void) transpose_X;
			gemv(y, 0., 1., transpose_A, A, x);
		}

		///gemv: y = (y * beta) + alpha * A * x
		static void gemv(Vector &y,
						 const Scalar beta,
						 const Scalar &alpha,
						 const bool transpose_A,
						 const Matrix &A,
						 const Vector &x);
		
		static void gemv(Vector &y,
						 const Scalar beta,
						 const Scalar &alpha,
						 const bool transpose_A,
						 const CRSMatrix<Scalar> &A,
						 const Vector &x);
		
		//[blas 3]

		inline static void multiply(Matrix &C,
						 const bool transpose_A,
						 const Matrix &A,
						 const bool transpose_B,
						 const Matrix &B)
		{
			gemm(C, 0., 1., transpose_A, A, transpose_B, B);
		}

		inline static void multiply(Matrix &C,
						 const bool transpose_A,
						 const Vector &A,
						 const bool transpose_B,
						 const Matrix &B)
		{
			gemm(C, 0., 1., transpose_A, A, transpose_B, B);
		}
		
		///gemm: C = beta * C + alpha * A * B
		static void gemm(Matrix &C,
						 const Scalar beta,
						 const Scalar &alpha,
						 const bool transpose_A,
						 const Matrix &A,
						 const bool transpose_B,
						 const Matrix &B);

        //if A is a Vector we still want to use this
       static void gemm(Matrix &C,
                         const Scalar beta,
                         const Scalar &alpha,
                         const bool transpose_A,
                         const Vector &A,
                         const bool transpose_B,
                         const Matrix &B);
		
		static void gemm(Matrix &C,
						 const Scalar beta,
						 const Scalar &alpha,
						 const bool transpose_A,
						 const CRSMatrix<Scalar> &A,
						 const bool transpose_B,
						 const CRSMatrix<Scalar> &B);
		
		//[unary]
		template<class Operation>
		static void apply_unary(Vector &result, const Operation &op, const Vector &vec) {
			if(vec.size() != result.size()) {
				result.resize(vec.size());
			}
			
			for(SizeType i = 0; i < vec.size(); ++i) {
				result[i] = op.template apply<Scalar>(vec[i]);
			}
		}
		
		static void apply_unary(Vector &result, const Minus &, const Vector &vec) {
			if(vec.size() != result.size()) {
				result.resize(vec.size());
			}
			
			for(SizeType i = 0; i < vec.size(); ++i) {
				result[i] = -vec[i];
			}
		}

		template<class Operation>
		static void apply_unary(Matrix &result, const Operation &op, const Matrix &mat) {
			if(result.rows() != mat.rows() || result.cols() != mat.cols()) {
				result.resize(mat.rows(), mat.cols());
			}

			for(SizeType i = 0; i < mat.size(); ++i) {
				result.entries()[i] = op.template apply<Scalar>(mat.entries()[i]);
			}
		}

		static void apply_unary(Matrix &result, const Minus &, const Matrix &mat) {
			if(result.rows() != mat.rows() || result.cols() != mat.cols()) {
				result.resize(mat.rows(), mat.cols());
			}

			for(SizeType i = 0; i < mat.size(); ++i) {
				result.entries()[i] = -mat.entries()[i];
			}
		}
		
		
		//[binary]
		template<class VectorT>
		inline static void apply_binary(Vector &result, const Scalar scaleFactor, const Multiplies &, VectorT &&v)
		{
			const int n = v.size();
			const int incx = 1;
			
			result = std::forward<VectorT>(v);
			assert(!result.empty());
			dscal_wrapper(&n, &scaleFactor, &result[0], &incx);
		}
		
		template<class MatrixT>
		inline static void apply_binary(Matrix &result, const Scalar scaleFactor, const Multiplies &, MatrixT &&m)
		{
			const int n = m.rows() * m.cols();
			const int incx = 1;
			result = std::forward<MatrixT>(m);
			dscal_wrapper(&n, &scaleFactor, result.ptr(), &incx);
		}

        template<class MatrixT>
        inline static void apply_binary(CRSMatrix<Scalar> &result, const Scalar scaleFactor, const Multiplies &, MatrixT &&m)
        {
            const int n = m.entries().size();
            const int incx = 1;
            result = std::forward<MatrixT>(m);
            dscal_wrapper(&n, &scaleFactor, ptr(result.entries()), &incx);
        }
		
		template<class Operation>
		inline static void apply_binary(Vector &result, const Vector &left, const Operation &op, const Vector &right) {
			if(result.size() != left.size()) {
				result.resize(left.size());
			}
			
			for (SizeType i = 0; i < left.size(); ++i) {
				result[i] = op.template apply<Scalar>(left[i], right[i]);
			}
		}
		

		static void apply_binary(Vector &result, const Reciprocal<Scalar> &reciprocal, const Vector &vec);
        static void apply_binary(Vector &result, const Vector &vec, const Multiplies &, const Scalar value);

		static void apply_binary(Vector &result, const Vector &vec,  const Divides &, const Scalar value);
		static void apply_binary(Vector &result, const Vector &left, const Minus &, const Vector &right);
		static void apply_binary(Vector &result, const Vector &left, const Plus &, const Vector &right);
		static void apply_binary(Vector &result, const Matrix &left, const Multiplies &, const Vector &right);
        static void apply_binary(Vector &result, const CRSMatrix<Scalar> &left, const Multiplies &, const Vector &right);

		static void apply_binary(Matrix &result, const Matrix &left, const Plus &, const Matrix &right);
		static void apply_binary(Matrix &result, const Matrix &left, const Minus &, const Matrix &right);
		static void apply_binary(Matrix &result, const Matrix &mat,  const Multiplies &, const Scalar &value);
		static void apply_binary(Matrix &result, const Matrix &mat,  const Divides &, const Scalar &value);
		static void apply_binary(Matrix &result, const Matrix &left, const Multiplies &, const Matrix &right);
		
		//[advanced]
		static void kronecker_product(Matrix &result, const Vector &left, const Vector &right);
		
		static void diag(Matrix &left, const Vector &right);
		static void diag(Vector &left, const Matrix &right);
		static void diag(Matrix &left, const Matrix &right);

		// static void diag(CRSMatrix<Scalar> &left, const Vector &right);
		static void diag(Vector &left, const CRSMatrix<Scalar>  &right);
		
		static void mat_diag_shift(Matrix &left, const Scalar diag_factor);
		static void diag_scale_right(Matrix &result, const Matrix &m, const Vector &diag)
		{
			result.resize(m.rows(), m.cols());
			ASSERT(diag.size() == m.cols() && "sizes are not compatible");
			
			for(SizeType i = 0; i < m.rows(); ++i) {
				for(SizeType j = 0; j < m.cols(); ++j) {
					result.set(i, j, m.get(i, j) * diag[j]);
				}
			}
		}
		
		static void diag_scale_left(Matrix &result, const Vector &diag, const Matrix &m)
		{
			result.resize(m.rows(), m.cols());
			
			const SizeType n = diag.size();
			ASSERT(n == m.rows() && "sizes are not compatible");
			
			if(n != m.rows()) {
				handle_error(-1);
				return;
			}
			
			for(SizeType i = 0; i < m.rows(); ++i) {
				for(SizeType j = 0; j < m.cols(); ++j) {
					result.set(i, j, m.get(i, j) * diag[i]);
				}
			}
		}
		//[selections]
		static void select(Vector &left,
						   const Vector &right,
						   const std::vector<SizeType> &index);
		
		static void select(Matrix &left,
						   const Matrix &right,
						   const std::vector<SizeType> &row_index,
						   const std::vector<SizeType> &col_index);
		
		
		static void assign_from_range(Matrix &left, const Matrix &right, const Range &row_range, const Range &col_range);
		static void assign_from_range(Vector &left, const Vector &right, const Range &row_range, const Range & /*col_range*/);
		static void assign_to_range(Matrix &left, const Matrix &right, const Range &row_range, const Range &col_range);
		static void assign_to_range(Matrix &left, const Identity &, const Range &row_range, const Range &col_range);
		static void assign_to_range(Vector &left, const Vector &right, const Range &row_range, const Range & /*col_range*/);
		
		//[misc]
		//FIXME what about this stuff??
		//------------------------------------------TODO-------------------------------- //
		template<class Tensor>
		static void monitor(const long &, Tensor &)
		{
			std::cout<<"BLASBackend:: monitor needs to be implemented  \n";
		}
		
		template<class Tensor>
		static void set_zero_rows(Tensor & /*Mat_A */, const std::vector<int> & /*index*/)
		{
			std::cout<<"BLASBackend:: set_zero_rows needs to be implemented  \n";
			
		}
		
		template<class MatTensor, class VecTensor>
		static void apply_BC_to_system(MatTensor & /*A*/, VecTensor& /*x*/, VecTensor& /*rhs*/, const std::vector<int> & /*index*/)
		{
			std::cout<<"BLASBackend:: apply_BC_to_system needs to be implemented  \n";
		}
		
		//-------------------------------------------------------------------------------------- //
		
		
		////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////
	};
	
	
	template<typename T>
	class Backend<T, BLAS> : public BLASBackend {
	public:
		inline static Backend &Instance()
		{
			static Backend instance;
			return instance;
		}
		
		BackendInfo &info()
		{
			return info_;
		}
		
	private:
		BackendInfo info_;
		
		Backend()
		{
			info_.set_name("blas");
		}
	};
}
#endif //utopia_utopia_BLASBACKEND_HPP
