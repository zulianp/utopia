#ifndef UTOPIA_TRILINOSBACKEND_HPP
#define UTOPIA_TRILINOSCBACKEND_HPP

#include "utopia_trilinos_Traits.hpp"

#include "utopia_Core.hpp"
#include "utopia_Factory.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_Base.hpp"
#include "utopia_ScalarBackend.hpp"
//#include "utopia_trilinos_Vector.hpp"
//#include "utopia_trilinos_Matrix.hpp"

#include <utility>

namespace utopia {	
	class TrilinosBackend : public ScalarBackend<double>  {
	public:
		typedef double Scalar;
		typedef TpetraVector Vector;
		typedef TpetraMatrix Matrix;

/*		using ScalarBackend<Scalar>::apply_binary;
		using ScalarBackend<Scalar>::axpy;
		
		template<class LorRValueMatrix>
		static void assign(TpetraMatrix &left, LorRValueMatrix &&right)
		{
			left = std::forward<LorRValueMatrix>(right);
		}
		
		template<class LorRValueVector>
		static void assign(TpetraVector &left, LorRValueVector &&right)
		{
			left = std::forward<LorRValueVector>(right);
		}

		static void assign_transposed(TpetraMatrix &left, const TpetraMatrix &right);
		
		static void clear(TpetraMatrix &mat);*/
				
/*		static void convert(Vec vec, TpetraVector &wrapper);
		static void convert(Mat mat, TpetraMatrix &wrapper);
		static void convert(Mat mat, PETScSparseMatrix &wrapper);
		static void convert(const TpetraVector &wrapper, Vec vec);
		static void convert(const TpetraMatrix &wrapper, Mat mat);
		static void convert(const PETScSparseMatrix &wrapper, Mat mat);
		static void wrap(Mat mat, PETScSparseMatrix &wrapper);
		static void wrap(Vec vec, TpetraVector &wrapper);
		
		static Range range(const TpetraVector &v);*/
		static Range row_range(const TpetraMatrix &m);
		static Range col_range(const TpetraMatrix &m);
		
		static void size(const TpetraMatrix &m, Size &size);
		static void size(const TpetraVector &m, Size &size);
		static void local_size(const TpetraVector &m, Size &size);
		static void local_size(const TpetraMatrix &m, Size &size);

		static void resize(TpetraMatrix &mat, const Size &size);
		static void resize(TpetraMatrix &mat, const Size &local_s, const Size &global_s);
		static void resize(TpetraVector &mat, const Size &size);
		static void resize(TpetraVector &vec, const Size &s_local, const Size &s_global);
		
		//[io]
		// read matrix
		static bool read(const std::string &path, TpetraMatrix &Mat_A);
		// write matrix
		static bool write(const std::string &path, const TpetraMatrix &Mat_A);
		
		// read vector
		static bool read(const std::string &path, TpetraVector &Vec_A);

		// write vector
		static bool write(const std::string &path, const TpetraVector &Vec_A);
		
		// monitor for cyrill 
		static void monitor(const long &it, TpetraMatrix &Mat_A); 
		static void monitor(const long &it, TpetraVector &Vec_A); 


		static Scalar get_global_nnz(TpetraMatrix &Mat_A); 
		static Scalar get_local_nnz(TpetraMatrix &Mat_A); 

		// call to MatZeroRows
		static void set_zero_rows(TpetraMatrix &Mat_A, const std::vector<int> &index); 
		static void apply_BC_to_system(TpetraMatrix & A, TpetraVector& x, TpetraVector& rhs, const std::vector<int> &index); 


		//[builders]
		template<class Tensor>
		static void build(Tensor &t, const Size &s, const Resize &)
		{
			build(t, s, Zeros());
		}
		
		static void build(TpetraMatrix &m, const Size &size, const Identity &);
		//static void build(PETScSparseMatrix &m, const Size &size, const Identity &);
		static void build(TpetraMatrix &m, const Size &size, const LocalIdentity &);
		/*static void build(PETScSparseMatrix &m, const Size &size, const LocalIdentity &);
		static void build(PETScSparseMatrix &m, const Size &size, const NNZ<int> &nnz);
		static void build(PETScSparseMatrix &m, const Size &size, const LocalNNZ<int> &nnz);
		static void build(PETScSparseMatrix &m, const Size &size, const LocalRowNNZ<int> &nnz);*/
		static void build(TpetraMatrix  &m, const Size &size, const LocalNNZ<int> & /*nnz */);
		static void build(TpetraMatrix  &m, const Size &size, const NNZ<int> &/*nnz*/);
		static void build(TpetraMatrix &m, const Size &size, const Zeros &);
		static void build(TpetraVector &v, const Size &size, const Zeros &);
		static void build(TpetraMatrix &m, const Size &size, const LocalZeros &);
		static void build(TpetraVector &v, const Size &size, const LocalZeros &);
		static void build(TpetraMatrix &m, const Size &size, const Values<Scalar> &values);
		static void build(TpetraVector &v, const Size &local_size, const Size &&global_size, const Values<Scalar> &values);
		static void build(TpetraVector &v, const Size &size, const Values<Scalar> &values);
		static void build(TpetraMatrix &m, const Size &size, const LocalValues<Scalar> &values);
		static void build(TpetraVector &v, const Size &size, const LocalValues<Scalar> &values);
		
		static void set(TpetraVector &v, const int index, Scalar value);
		static void add(TpetraVector &v, const int index, Scalar value);
		
		static void set(TpetraVector &v, const std::vector<int> indices, const std::vector<Scalar> values);
		static void set(TpetraMatrix &v, const int row, const int col, Scalar value);
		//static void set(PETScSparseMatrix &v, const int row, const int col, Scalar value);
		static void add(TpetraMatrix &v, const int row, const int col, Scalar value);
		

		//[host/device locks]
		template<class Tensor>
		static void read_lock(const Tensor &) {}
		
		template<class Tensor>
		static void read_unlock(const Tensor &) {}
		
		static void write_lock(TpetraVector &vec);
		static void write_unlock(TpetraVector &vec);
		static void write_lock(const TpetraMatrix &mat);
		static void write_unlock(const TpetraMatrix &mat);
		
		/*static void write_lock(const PETScSparseMatrix &mat);
		static void write_unlock(const PETScSparseMatrix &mat);*/
		
		static void set(TpetraMatrix &v, const std::vector<int> rows, const std::vector<int> cols, const std::vector<Scalar> values);
		static Scalar get(const TpetraVector &v, const int index);
		static Scalar get(const TpetraMatrix &v, const int row, const int col);

		//[unary]
		static void apply_unary(Vector &result, const Abs &, const Vector &vec);

		template<class Operation>
		static void apply_unary(Vector &result, const Operation &op, const Vector &v)
		{
			//Range r = range(v); //TODO

			Size gs, ls;
			size(v, gs);
			local_size(v, ls);

			/*VecType type;
			VecGetType(v.implementation(), &type);
			VecSetType(result.implementation(), type);
			VecSetSizes(result.implementation(), ls.get(0), gs.get(0));*/

			write_lock(result);
			read_lock(v);

			/*for(int i = r.begin(); i < r.end(); ++i) {
				Scalar value = get(v, i);
				VecSetValue(result.implementation(), i, op.template apply<Scalar>(value), INSERT_VALUES);
			}*/

			write_unlock(result);
			read_unlock(v);
		}

		static void apply_unary(Vector &result, const Vector &v, const Minus &)
		{
			//Range r = range(v);

			Size gs, ls;
			size(v, gs);
			local_size(v, ls);

			/*VecType type;
			VecGetType(v.implementation(), &type);
			VecSetType(result.implementation(), type);
			VecSetSizes(result.implementation(), ls.get(0), gs.get(0));*/

			write_lock(result);
			read_lock(v);

			/*for(int i = r.begin(); i < r.end(); ++i) {
				Scalar value = get(v, i);
				VecSetValue(result.implementation(), i, -value, INSERT_VALUES);
			}*/

			write_unlock(result);
			read_unlock(v);
		}

		//[binary]
		static void apply_binary(TpetraVector &result, const Reciprocal<Scalar> &reciprocal, const TpetraVector &vec);
		static void apply_binary(Vector &result, const Scalar factor, const Multiplies &, const Vector &vec);
		static void apply_binary(Matrix &result, const Scalar factor, const Multiplies &, const Matrix &mat);
		static void apply_binary(TpetraVector &result, const TpetraMatrix &left, const Multiplies &, const TpetraVector &right);
		static void apply_binary(TpetraMatrix &result, const TpetraMatrix &left, const Multiplies &, const TpetraMatrix &right);
		static void apply_binary(TpetraVector &result, const TpetraVector &left, const EMultiplies &, const TpetraVector &right);
		static void apply_binary(TpetraVector &result, const TpetraVector &left, const Divides &, const TpetraVector &right);
		static void apply_binary(TpetraVector &result, const TpetraVector &left, const Min &, const TpetraVector &right);
		static void apply_binary(TpetraVector &result, const TpetraVector &left, const Max &, const TpetraVector &right);

		template<class LeftTensor, class RightTensor, class ResultTensor>
		static void apply_binary(ResultTensor &result, LeftTensor &&left, const Plus &, RightTensor &&right) {
			if(&result == &left) {
				axpy(result, 1., right);
				return;
			}

			if(&result == &right) {
				axpy(result, 1., left);
				return;
			}

			result = std::forward<RightTensor>(right);
			axpy(result, 1., left);
		}
		
		template<class LeftTensor, class ResultTensor, class RightTensor>
		static void apply_binary(ResultTensor &result, LeftTensor &&left, const Minus &, RightTensor &&right) {
			assert(&result != &right);

			result = std::forward<LeftTensor>(left);
			axpy(result, -1., right);
		}

		static void allocate_apply_vec(TpetraVector &result, const TpetraVector &left, const TpetraVector &right);

		template<class Operation>
		static void apply_binary_generic(TpetraVector &result, const TpetraVector &left, const Operation &op, const TpetraVector &right)
		{
			allocate_apply_vec(result, left, right);

			const vector_type &l ;//= left.implementation();
			const vector_type &r ;//= right.implementation();
			vector_type &out ;//= result.implementation();

			/*auto ll = range(left);  TODO
			auto rr = range(right);*/

			/*assert(ll.extent() == rr.extent());
			assert(ll.begin() == rr.begin());*/

			read_lock(left);
			read_lock(right);
			write_lock(result);

			//for(int i = ll.begin(); i < ll.end(); ++i) {//TODO
				Scalar lv, rv;
				//VecGetValues(l, 1, &i, &lv);
				//VecGetValues(r, 1, &i, &rv);
				Scalar val = op.template apply<Scalar>(lv, rv);
			//	VecSetValues(out, 1, &i, &val, INSERT_VALUES); //TODO
			//}

			read_unlock(left);
			read_unlock(right);
			write_unlock(result);
		}

		//[specialized]
		static void mat_mult_add(TpetraVector &result, const TpetraMatrix &m, const TpetraVector &right, const TpetraVector &left);
		static void mat_mult_t_add(TpetraVector &result, const TpetraMatrix &m, const TpetraVector &right, const TpetraVector &left);		
		static void triple_product_ptap(TpetraMatrix &result, const TpetraMatrix &, const TpetraMatrix &); 
		static void triple_product(TpetraMatrix &result, const TpetraMatrix &, const TpetraMatrix &, const TpetraMatrix &); 


		static Scalar norm2(const TpetraVector &v);
		static Scalar norm2(const Matrix &m);
		static Scalar norm1(const TpetraVector &v);
		static Scalar norm_infty(const TpetraVector &v);
		static Scalar reduce(const TpetraVector &vec, const Plus &);
		static Scalar reduce(const TpetraMatrix &mat, const Plus &op);
		static Scalar reduce(const TpetraVector &, const Min &);
		static Scalar reduce(const TpetraMatrix &, const Min &);
		static Scalar reduce(const TpetraVector &, const Max &);
		static Scalar reduce(const TpetraMatrix &, const Max &);
		static Scalar dot(const Vector &left, const Vector &right);
		
		// get diagonal of matrix as vector
		static void diag(TpetraVector &vec, const TpetraMatrix &mat);
		//static void diag(PETScSparseMatrix &mat, const TpetraVector &vec);
		static void diag(TpetraMatrix &mat, const TpetraVector &vec);
		static void diag(TpetraMatrix &out, const TpetraMatrix &in);
        static void mat_diag_shift(TpetraMatrix &left, const Scalar diag_factor);
        

		static bool compare(const Vector &left, const Vector &right, const ApproxEqual &comp);
		static bool compare(const Matrix &left, const Matrix &right, const ApproxEqual &comp);
		
		
		static void kronecker_product(Matrix &result, const Vector &left, const Vector &right);
		
		//[selection]
		static void assign_from_range(TpetraVector &left, const TpetraVector &right, const Range &globalRowRange,
							 const Range &global_col_range);
		
		static void assign_to_range(TpetraMatrix &left, const TpetraMatrix &right, const Range &global_row_range,
						   const Range &global_col_range);
		
		static void assign_to_range( TpetraMatrix &left, const Identity &, const Range &global_row_range,
						   const Range &global_col_range);
		
		static void assign_to_range( TpetraVector &left, const TpetraVector &right, const Range &global_row_range,
						   const Range &global_col_range);

		static void assign_from_range(
			TpetraMatrix &left,
			const TpetraMatrix &right,
			const Range &globalRowRange,
			const Range &global_col_range);

	

		static void select(
					TpetraVector &left,
					const TpetraVector &right,
		      		const std::vector<int> &index);

		static void select(
					TpetraMatrix &left,
					const TpetraMatrix &right,
		      		const std::vector<int> &row_index,
		      		const std::vector<int> &col_index);


	
		// this is not very correct yet ...
		//static void build_local_diag_block(PETScSerialSparseMatrix &left, const PETScSparseMatrix &right);
		static void build_local_redistribute(TpetraVector &result, const TpetraVector &, const TpetraVector &); 

		static bool is_nan_or_inf(const TpetraVector &v); 
		static bool is_nan_or_inf(const TpetraMatrix &m); 


		template<class Tensor>
		static void read_and_write_lock(Tensor &t) {
			write_lock(t);
		}
		
		template<class Tensor>
		static void read_and_write_unlock(Tensor &t){
			write_unlock(t);
		}

		static void diag_scale_right(Matrix &result, const Matrix &m,    const Vector &diag);
		static void diag_scale_left(Matrix &result,  const Vector &diag, const Matrix &m);
		static void diag_scale_left(Vector &result,  const Vector &diag, const Vector &m);

		static Scalar trace(const Matrix &mat);
		static void apply_tensor_reduce(Vector &result, const Matrix &mat, const Plus &, const int dim);
		static void apply_tensor_reduce(Vector &result, const Matrix &mat, const Min &,  const int dim);
		static void apply_tensor_reduce(Vector &result, const Matrix &mat, const Max &,  const int dim);
		static void inverse(Matrix &result, const Matrix &mat);

		static void axpy(TpetraVector &y, const Scalar alpha, const TpetraVector &x);
		static void axpy(TpetraMatrix &y, const Scalar alpha, const TpetraMatrix &x);
		static void axpby(Vector &y, const Scalar alpha, const Vector &x, const Scalar &beta);

		inline static void multiply(
			TpetraMatrix &result,
			bool transpose_left,
			const TpetraMatrix &left,
			bool transpose_right,
			const TpetraMatrix &right)
		{
			gemm(result, 0.0, 1., transpose_left, left, transpose_right, right);
		}

		inline static void multiply(
			Vector &result,
			bool transpose_left,
			const TpetraMatrix &left,
			bool transpose_right,
			const Vector &right)
		{
			assert(!transpose_right); 
			(void) transpose_right;
			
			gemv(result, 0.0, 1., transpose_left, left, right);
		}

	private:

		static void gemm(
			TpetraMatrix &result,
			const Scalar beta,
			const Scalar alpha,
			bool transpose_left,
			const TpetraMatrix &left,
			bool transpose_right,
			const TpetraMatrix &right);

		static void gemv(
			Vector &y, 
			const Scalar beta,
			const Scalar alpha,
			bool transpose_A,
			const TpetraMatrix &A, 
			const Vector &x);
				
		static void select_aux(TpetraMatrix &left,
					const TpetraMatrix &right,
		      		const std::vector<int> &row_index,
		      		const std::vector<int> &col_index);

		static void par_assign_from_local_is(
			const std::vector<int> &remote_rows,
			const std::vector<int> &remote_cols,
			const int global_col_offset,
			const Range &local_col_range,
			const TpetraMatrix &right,
			TpetraMatrix &result);

		static void par_assign_from_local_range(
							const Range &local_row_range,
							const Range &local_col_range,
							const Range &global_col_range,
							const TpetraMatrix &right,
							TpetraMatrix &result);

		static void par_assign_from_global_range(
		const Range &global_row_range,
		const Range &global_col_range,
		const TpetraMatrix &right,
		TpetraMatrix &result);


		//unused
		static void vec_to_mat(Matrix &m, const Vector &v, const bool transpose);
	};
	
	template<>
	class Backend<double,TRILINOS> : public TrilinosBackend {  //TODO Trilinos second template parameter
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
			info_.set_name("trilinos");
		}
	};
}

#endif //UTOPIA_TRILINOSBACKEND_HPP

