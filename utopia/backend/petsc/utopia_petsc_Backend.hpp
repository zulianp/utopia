#ifndef UTOPIA_UTOPIA_PETSCBACKEND_HPP
#define UTOPIA_UTOPIA_PETSCBACKEND_HPP

#include "utopia_petsc_Traits.hpp"

#include "utopia_Core.hpp"
#include "utopia_Factory.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_Base.hpp"
#include "utopia_ScalarBackend.hpp"

#include <utility>

namespace utopia {	
	class PetscBackend : public ScalarBackend<PetscScalar>  {
	public:
		typedef PetscScalar Scalar;
		typedef PETScVector Vector;
		typedef PETScMatrix Matrix;

		class PetscArgs {
		public: 
			MPI_Comm comm;
			std::string name;

			template<class... Args>
			void parse(const Optional<Args...> &args)
			{
				args.each(*this);
			}

			template<class Any>
			inline static void parse_arg(const Any &) {}

			void parse_arg(const MPI_Comm &comm)
			{
				this->comm = comm;
			}

			void parse_arg(const std::string &name)
			{
				//std::cout << "name: " << name << std::endl;
				this->name = name;
			}

			PetscArgs()
			: comm(PetscBackend::default_communicator()), name()
			{}
		};

		template<class... Args>
		PetscArgs parse_args(const Optional<Args...> &opts)
		{
			PetscArgs args;
			args.parse(opts);
			return args;
		}
		
		using ScalarBackend<Scalar>::apply_binary;
		using ScalarBackend<Scalar>::axpy;
		
		template<class LorRValueMatrix>
		static void assign(PETScMatrix &left, LorRValueMatrix &&right)
		{
			left = std::forward<LorRValueMatrix>(right);
		}
		
		template<class LorRValueVector>
		static void assign(PETScVector &left, LorRValueVector &&right)
		{
			left = std::forward<LorRValueVector>(right);
		}

		static void assign_transposed(PETScMatrix &left, const PETScMatrix &right);
		
		static void clear(PETScMatrix &mat);
				
		static void convert(Vec vec, PETScVector &wrapper);
		static void convert(Mat mat, PETScMatrix &wrapper);
		static void convert(Mat mat, PETScSparseMatrix &wrapper);
		static void convert(const PETScVector &wrapper, Vec vec);
		static void convert(const PETScMatrix &wrapper, Mat mat);
		static void convert(const PETScSparseMatrix &wrapper, Mat mat);
		static void wrap(Mat mat, PETScSparseMatrix &wrapper);
		static void wrap(Vec vec, PETScVector &wrapper);
		
		static Range range(const PETScVector &v);
		static Range row_range(const PETScMatrix &m);
		static Range col_range(const PETScMatrix &m);
		
		static void size(const PETScMatrix &m, Size &size);
		static void size(const PETScVector &m, Size &size);
		static void local_size(const PETScVector &m, Size &size);
		static void local_size(const PETScMatrix &m, Size &size);

		static void resize(PETScMatrix &mat, const Size &size);
		static void resize(PETScMatrix &mat, const Size &local_s, const Size &global_s);
		static void resize(PETScVector &mat, const Size &size);
		static void resize(PETScVector &vec, const Size &s_local, const Size &s_global);
		
		//[io]
		// read matrix
		static bool read(const std::string &path, PETScMatrix &Mat_A, const PetscArgs &args = PetscArgs());
		// write matrix
		static bool write(const std::string &path, const PETScMatrix &Mat_A);
		
		// read vector
		static bool read(const std::string &path, PETScVector &Vec_A, const PetscArgs &args = PetscArgs());

		// write vector
		static bool write(const std::string &path, const PETScVector &Vec_A);
		
		// monitoring functions for iterative solvers (Cyrill)
		static void monitor(const long &iteration, PETScMatrix &m); 
		static void monitor(const long &iteration, PETScVector &v); 

		static Scalar get_global_nnz(PETScMatrix &Mat_A); 
		static Scalar get_local_nnz(PETScMatrix &Mat_A); 

		// call to MatZeroRows
		static void set_zero_rows(PETScMatrix &Mat_A, const std::vector<int> &index); 
		static void apply_BC_to_system(PETScMatrix & A, PETScVector& x, PETScVector& rhs, const std::vector<int> &index); 


		//[builders]
		template<class Tensor>
		static void build(Tensor &t, const Size &s, const Resize &, const PetscArgs &opts = PetscArgs())
		{
			build(t, s, Zeros(), opts);
		}
		
		static void build(PETScMatrix &m, const Size &size, const Identity &, const PetscArgs &opts = PetscArgs());
		static void build(PETScSparseMatrix &m, const Size &size, const Identity &, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts = PetscArgs());
		static void build(PETScSparseMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts = PetscArgs());
		static void build(PETScSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz, const PetscArgs &opts = PetscArgs());
		static void build(PETScSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix  &m, const Size &size, const LocalNNZ<PetscInt> & /*nnz */, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix  &m, const Size &size, const NNZ<PetscInt> &/*nnz*/, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix &m, const Size &size, const Zeros &, const PetscArgs &opts = PetscArgs());
		static void build(PETScVector &v, const Size &size, const Zeros &, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix &m, const Size &size, const LocalZeros &, const PetscArgs &opts = PetscArgs());
		static void build(PETScVector &v, const Size &size, const LocalZeros &, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix &m, const Size &size, const Values<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PETScVector &v, const Size &local_size, const Size &&global_size, const Values<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PETScVector &v, const Size &size, const Values<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PETScMatrix &m, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PETScVector &v, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts = PetscArgs());
		
		static void set(PETScVector &v, const PetscInt index, Scalar value);
		static void add(PETScVector &v, const PetscInt index, Scalar value);
		
		static void set(PETScVector &v, const std::vector<PetscInt> &indices, const std::vector<Scalar> &values);
		static void set(PETScMatrix &v, const PetscInt row, const PetscInt col, Scalar value);
		static void set(PETScSparseMatrix &v, const PetscInt row, const PetscInt col, Scalar value);
		static void add(PETScMatrix &m, const PetscInt row, const PetscInt col, Scalar value);
		

		//[host/device locks]
		template<class Tensor>
		static void read_lock(const Tensor &) {}
		
		template<class Tensor>
		static void read_unlock(const Tensor &) {}
		
		static void write_lock(PETScVector &vec);
		static void write_unlock(PETScVector &vec);
		static void write_lock(const PETScMatrix &mat);
		static void write_unlock(const PETScMatrix &mat);
		
		static void write_lock(const PETScSparseMatrix &mat);
		static void write_unlock(const PETScSparseMatrix &mat);
		
		static void set(PETScMatrix &v, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);
		
		template<typename T>
		static void add_matrix(PETScMatrix &m, const std::vector<T> &rows, const std::vector<T> &cols, const std::vector<Scalar> &values)
		{
			std::vector<PetscInt> petsc_rows, petsc_cols;
			petsc_rows.insert(petsc_rows.end(), rows.begin(), rows.end());
			petsc_cols.insert(petsc_cols.end(), cols.begin(), cols.end());
			PetscBackend::add_matrix(m, petsc_rows, petsc_cols, values);
		}


		static void add_matrix(PETScMatrix &v, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);
		static void add_matrix(PETScSparseMatrix &v, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);

		static Scalar get(const PETScVector &v, const PetscInt index);
		static Scalar get(const PETScMatrix &v, const PetscInt row, const PetscInt col);

		//[unary]
		static void apply_unary(Vector &result, const Abs &, const Vector &vec);

		template<class Operation>
		static void apply_unary(Vector &result, const Operation &op, const Vector &v)
		{
			Range r = range(v);

			Size gs, ls;
			size(v, gs);
			local_size(v, ls);

			if(result.is_null()) {
				VecCreate(v.communicator(), &result.implementation());
			} else if(result.communicator() != v.communicator()) {
				VecDestroy(&result.implementation());
				VecCreate(v.communicator(), &result.implementation());
			}

			VecType type;
			VecGetType(v.implementation(), &type);
			VecSetType(result.implementation(), type);
			VecSetSizes(result.implementation(), ls.get(0), gs.get(0));

			write_lock(result);
			read_lock(v);

			for(PetscInt i = r.begin(); i < r.end(); ++i) {
				Scalar value = get(v, i);
				VecSetValue(result.implementation(), i, op.template apply<Scalar>(value), INSERT_VALUES);
			}

			write_unlock(result);
			read_unlock(v);
		}

		static void apply_unary(Vector &result, const Vector &v, const Minus &)
		{
			Range r = range(v);

			Size gs, ls;
			size(v, gs);
			local_size(v, ls);

			if(result.is_null()) {
				VecCreate(v.communicator(), &result.implementation());
			} else if(result.communicator() != v.communicator()) {
				VecDestroy(&result.implementation());
				VecCreate(v.communicator(), &result.implementation());
			}

			VecType type;
			VecGetType(v.implementation(), &type);
			VecSetType(result.implementation(), type);
			VecSetSizes(result.implementation(), ls.get(0), gs.get(0));

			write_lock(result);
			read_lock(v);

			for(PetscInt i = r.begin(); i < r.end(); ++i) {
				Scalar value = get(v, i);
				VecSetValue(result.implementation(), i, -value, INSERT_VALUES);
			}

			write_unlock(result);
			read_unlock(v);
		}

		//[binary]
		static void apply_binary(PETScVector &result, const Reciprocal<Scalar> &reciprocal, const PETScVector &vec);
		static void apply_binary(Vector &result, const Scalar factor, const Multiplies &, const Vector &vec);
		static void apply_binary(Matrix &result, const Scalar factor, const Multiplies &, const Matrix &mat);
		static void apply_binary(PETScVector &result, const PETScMatrix &left, const Multiplies &, const PETScVector &right);
		static void apply_binary(PETScMatrix &result, const PETScMatrix &left, const Multiplies &, const PETScMatrix &right);
		static void apply_binary(PETScVector &result, const PETScVector &left, const EMultiplies &, const PETScVector &right);
		static void apply_binary(PETScVector &result, const PETScVector &left, const Divides &, const PETScVector &right);
		static void apply_binary(PETScVector &result, const PETScVector &left, const Min &, const PETScVector &right);
		static void apply_binary(PETScVector &result, const PETScVector &left, const Max &, const PETScVector &right);

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

		static void allocate_apply_vec(PETScVector &result, const PETScVector &left, const PETScVector &right);

		template<class Operation>
		static void apply_binary_generic(PETScVector &result, const PETScVector &left, const Operation &op, const PETScVector &right)
		{
			allocate_apply_vec(result, left, right);

			const Vec &l = left.implementation();
			const Vec &r = right.implementation();
			Vec &out = result.implementation();

			auto ll = range(left);

			assert(ll.extent() == range(right).extent());
			assert(ll.begin()  == range(right).begin());

			read_lock(left);
			read_lock(right);
			write_lock(result);

			for(PetscInt i = ll.begin(); i < ll.end(); ++i) {
				Scalar lv, rv;
				VecGetValues(l, 1, &i, &lv);
				VecGetValues(r, 1, &i, &rv);
				Scalar val = op.template apply<Scalar>(lv, rv);
				VecSetValues(out, 1, &i, &val, INSERT_VALUES);
			}

			read_unlock(left);
			read_unlock(right);
			write_unlock(result);
		}

		//[specialized]
		static void mat_mult_add(PETScVector &result, const PETScMatrix &m, const PETScVector &right, const PETScVector &left);
		static void mat_mult_t_add(PETScVector &result, const PETScMatrix &m, const PETScVector &right, const PETScVector &left);		
		static void triple_product_ptap(PETScMatrix &result, const PETScMatrix &, const PETScMatrix &); 
		static void triple_product(PETScMatrix &result, const PETScMatrix &, const PETScMatrix &, const PETScMatrix &); 


		static Scalar norm2(const PETScVector &v);
		static Scalar norm2(const Matrix &m);
		static Scalar norm1(const PETScVector &v);
		static Scalar norm_infty(const PETScVector &v);
		static Scalar reduce(const PETScVector &vec, const Plus &);
		static Scalar reduce(const PETScMatrix &mat, const Plus &op);
		static Scalar reduce(const PETScVector &, const Min &);
		static Scalar reduce(const PETScMatrix &, const Min &);
		static Scalar reduce(const PETScVector &, const Max &);
		static Scalar reduce(const PETScMatrix &, const Max &);
		static Scalar dot(const Vector &left, const Vector &right);
		
		// get diagonal of matrix as vector
		static void diag(PETScVector &vec, const PETScMatrix &mat);
		static void diag(PETScSparseMatrix &mat, const PETScVector &vec);
		static void diag(PETScMatrix &mat, const PETScVector &vec);
		static void diag(PETScMatrix &out, const PETScMatrix &in);
        static void mat_diag_shift(PETScMatrix &left, const Scalar diag_factor);
        

		static bool compare(const Vector &left, const Vector &right, const ApproxEqual &comp);
		static bool compare(const Matrix &left, const Matrix &right, const ApproxEqual &comp);
		
		
		static void kronecker_product(Matrix &result, const Vector &left, const Vector &right);
		
		//[selection]
		static void assign_from_range(PETScVector &left, const PETScVector &right, const Range &globalRowRange,
							 const Range &global_col_range);
		
		static void assign_to_range(PETScMatrix &left, const PETScMatrix &right, const Range &global_row_range,
						   const Range &global_col_range);
		
		static void assign_to_range( PETScMatrix &left, const Identity &, const Range &global_row_range,
						   const Range &global_col_range);
		
		static void assign_to_range( PETScVector &left, const PETScVector &right, const Range &global_row_range,
						   const Range &global_col_range);

		static void assign_from_range(
			PETScMatrix &left,
			const PETScMatrix &right,
			const Range &globalRowRange,
			const Range &global_col_range);

	

		static void select(
					PETScVector &left,
					const PETScVector &right,
		      		const std::vector<PetscInt> &index);

		static void select(
					PETScMatrix &left,
					const PETScMatrix &right,
		      		const std::vector<PetscInt> &row_index,
		      		const std::vector<PetscInt> &col_index);


	
		// this is not very correct yet ...
		static void build_local_diag_block(PETScSerialSparseMatrix &left, const PETScSparseMatrix &right);
		static void build_local_redistribute(PETScVector &result, const PETScVector &, const PETScVector &); 

		static bool is_nan_or_inf(const PETScVector &v); 
		static bool is_nan_or_inf(const PETScMatrix &m); 


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

		static void axpy(PETScVector &y, const Scalar alpha, const PETScVector &x);
		static void axpy(PETScMatrix &y, const Scalar alpha, const PETScMatrix &x);
		static void axpby(Vector &y, const Scalar alpha, const Vector &x, const Scalar &beta);

		inline static void multiply(
			PETScMatrix &result,
			bool transpose_left,
			const PETScMatrix &left,
			bool transpose_right,
			const PETScMatrix &right)
		{
			gemm(result, 0.0, 1., transpose_left, left, transpose_right, right);
		}

		inline static void multiply(
			Vector &result,
			bool transpose_left,
			const PETScMatrix &left,
			bool transpose_right,
			const Vector &right)
		{
			assert(!transpose_right); 
			(void) transpose_right;
			
			gemv(result, 0.0, 1., transpose_left, left, right);
		}

		static void scale(Vector &result, const Scalar scale_factor);
		static void scale(Matrix &result, const Scalar scale_factor);

		static void vec_create_parallel(MPI_Comm comm, PetscInt n_local, PetscInt n_global, Vec *vec);
		static void vec_repurpose(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global, Vec *vec);

		static void sparse_mat_create_parallel(
				MPI_Comm comm,
				PetscInt rows_local,
				PetscInt cols_local,
				PetscInt rows_global,
				PetscInt cols_global,
				PetscInt d_nnz,
				PetscInt o_nnz,
				Mat *mat);


		static void dense_mat_create_parallel(
				MPI_Comm comm,
				PetscInt rows_local,
				PetscInt cols_local,
				PetscInt rows_global,
				PetscInt cols_global,
				Mat *mat);

		inline static bool check_error(const PetscInt err) 
		{
			return PETScError::Check(err);
		}

		inline static MatType parallel_sparse_matrix_type()
		{
			return MATAIJ;
		}

		inline static MatType parallel_dense_matrix_type()
		{
			return MATDENSE;
		}

		inline static VecType parallel_vector_type()
		{
			return VECSTANDARD;
		}

	private:

		static void gemm(
			PETScMatrix &result,
			const Scalar beta,
			const Scalar alpha,
			bool transpose_left,
			const PETScMatrix &left,
			bool transpose_right,
			const PETScMatrix &right);

		static void gemv(
			Vector &y, 
			const Scalar beta,
			const Scalar alpha,
			bool transpose_A,
			const PETScMatrix &A, 
			const Vector &x);
				
		static void select_aux(PETScMatrix &left,
					const PETScMatrix &right,
		      		const std::vector<PetscInt> &row_index,
		      		const std::vector<PetscInt> &col_index);

		static void par_assign_from_local_is(
			const std::vector<PetscInt> &remote_rows,
			const std::vector<PetscInt> &remote_cols,
			const PetscInt global_col_offset,
			const Range &local_col_range,
			const PETScMatrix &right,
			PETScMatrix &result);

		static void par_assign_from_local_range(
							const Range &local_row_range,
							const Range &local_col_range,
							const Range &global_col_range,
							const PETScMatrix &right,
							PETScMatrix &result);

		static void par_assign_from_global_range(
		const Range &global_row_range,
		const Range &global_col_range,
		const PETScMatrix &right,
		PETScMatrix &result);


		//unused
		static void vec_to_mat(Matrix &m, const Vector &v, const bool transpose);

		inline static MPI_Comm default_communicator() {
			return PETSC_COMM_WORLD;
		}

		static void apply_args(const PetscArgs &args, Matrix &m);
		static void apply_args(const PetscArgs &args, Vector &m);
	};
	
	template<>
	class Backend<PetscScalar, PETSC> : public PetscBackend {
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
			info_.set_name("petsc");
		}
	};
}

#endif //UTOPIA_UTOPIA_PETSCBACKEND_HPP

