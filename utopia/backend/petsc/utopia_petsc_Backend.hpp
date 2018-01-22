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
		typedef PetscVector Vector;
		typedef PetscMatrix Matrix;

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
		static void assign(PetscMatrix &left, LorRValueMatrix &&right)
		{
			left = std::forward<LorRValueMatrix>(right);
		}
		
		template<class LorRValueVector>
		static void assign(PetscVector &left, LorRValueVector &&right)
		{
			left = std::forward<LorRValueVector>(right);
		}

		static void assign_transposed(PetscMatrix &left, const PetscMatrix &right);
		
		static void clear(PetscMatrix &mat);
				
		static void convert(Vec vec, PetscVector &wrapper);
		static void convert(Mat mat, PetscMatrix &wrapper);
		static void convert(Mat mat, PetscSparseMatrix &wrapper);
		static void convert(const PetscVector &wrapper, Vec vec);
		static void convert(const PetscMatrix &wrapper, Mat mat);
		static void convert(const PetscSparseMatrix &wrapper, Mat mat);
		static void wrap(Mat mat, PetscSparseMatrix &wrapper);
		static void wrap(Vec vec, PetscVector &wrapper);
		
		static Range range(const PetscVector &v);
		static Range row_range(const PetscMatrix &m);
		static Range col_range(const PetscMatrix &m);
		
		static void size(const PetscMatrix &m, Size &size);
		static void size(const PetscVector &m, Size &size);
		static void local_size(const PetscVector &m, Size &size);
		static void local_size(const PetscMatrix &m, Size &size);

		static void resize(PetscMatrix &mat, const Size &size);
		static void resize(PetscMatrix &mat, const Size &local_s, const Size &global_s);
		static void resize(PetscVector &mat, const Size &size);
		static void resize(PetscVector &vec, const Size &s_local, const Size &s_global);
		
		//[io]
		// read matrix
		static bool read(const std::string &path, PetscMatrix &Mat_A, const PetscArgs &args = PetscArgs());
		// write matrix
		static bool write(const std::string &path, const PetscMatrix &Mat_A);
		
		// read vector
		static bool read(const std::string &path, PetscVector &Vec_A, const PetscArgs &args = PetscArgs());

		// write vector
		static bool write(const std::string &path, const PetscVector &Vec_A);
		
		// monitoring functions for iterative solvers (Cyrill)
		static void monitor(const long &iteration, PetscMatrix &m); 
		static void monitor(const long &iteration, PetscVector &v); 

		static Scalar get_global_nnz(PetscMatrix &Mat_A); 
		static Scalar get_local_nnz(PetscMatrix &Mat_A); 

		// call to MatZeroRows
		static void set_zero_rows(PetscMatrix &Mat_A, const std::vector<int> &index); 
		static void apply_BC_to_system(PetscMatrix & A, PetscVector& x, PetscVector& rhs, const std::vector<int> &index); 


		//[builders]
		template<class Tensor>
		static void build(Tensor &t, const Size &s, const Resize &, const PetscArgs &opts = PetscArgs())
		{
			build(t, s, Zeros(), opts);
		}
		
		static void build(PetscMatrix &m, const Size &size, const Identity &, const PetscArgs &opts = PetscArgs());
		static void build(PetscSparseMatrix &m, const Size &size, const Identity &, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts = PetscArgs());
		static void build(PetscSparseMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts = PetscArgs());
		static void build(PetscSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz, const PetscArgs &opts = PetscArgs());
		static void build(PetscSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix  &m, const Size &size, const LocalNNZ<PetscInt> & /*nnz */, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix  &m, const Size &size, const NNZ<PetscInt> &/*nnz*/, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix &m, const Size &size, const Zeros &, const PetscArgs &opts = PetscArgs());
		static void build(PetscVector &v, const Size &size, const Zeros &, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix &m, const Size &size, const LocalZeros &, const PetscArgs &opts = PetscArgs());
		static void build(PetscVector &v, const Size &size, const LocalZeros &, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix &m, const Size &size, const Values<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PetscVector &v, const Size &local_size, const Size &&global_size, const Values<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PetscVector &v, const Size &size, const Values<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PetscMatrix &m, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts = PetscArgs());
		static void build(PetscVector &v, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts = PetscArgs());
		
		static void set(PetscVector &v, const PetscInt index, Scalar value);
		static void add(PetscVector &v, const PetscInt index, Scalar value);
		
		static void set(PetscVector &v, const std::vector<PetscInt> &indices, const std::vector<Scalar> &values);
		static void set(PetscMatrix &v, const PetscInt row, const PetscInt col, Scalar value);
		static void set(PetscSparseMatrix &v, const PetscInt row, const PetscInt col, Scalar value);
		static void add(PetscMatrix &m, const PetscInt row, const PetscInt col, Scalar value);
		

		//[host/device locks]
		template<class Tensor>
		static void read_lock(const Tensor &) {}
		
		template<class Tensor>
		static void read_unlock(const Tensor &) {}
		
		static void write_lock(PetscVector &vec);
		static void write_unlock(PetscVector &vec);
		static void write_lock(PetscMatrix &mat);
		static void write_unlock(PetscMatrix &mat);
		
		static void set(PetscMatrix &v, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);
		
		template<typename T>
		static void add_matrix(PetscMatrix &m, const std::vector<T> &rows, const std::vector<T> &cols, const std::vector<Scalar> &values)
		{
			std::vector<PetscInt> petsc_rows, petsc_cols;
			petsc_rows.insert(petsc_rows.end(), rows.begin(), rows.end());
			petsc_cols.insert(petsc_cols.end(), cols.begin(), cols.end());
			PetscBackend::add_matrix(m, petsc_rows, petsc_cols, values);
		}

		template<typename T>
		static void set_matrix(PetscMatrix &m, const std::vector<T> &rows, const std::vector<T> &cols, const std::vector<Scalar> &values)
		{
			std::vector<PetscInt> petsc_rows, petsc_cols;
			petsc_rows.insert(petsc_rows.end(), rows.begin(), rows.end());
			petsc_cols.insert(petsc_cols.end(), cols.begin(), cols.end());
			PetscBackend::set_matrix(m, petsc_rows, petsc_cols, values);
		}

		static void add_matrix(PetscMatrix &m, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);
		static void add_matrix(PetscSparseMatrix &m, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);

		static void set_matrix(PetscMatrix &m, const std::vector<PetscInt> &rows, const std::vector<PetscInt> &cols, const std::vector<Scalar> &values);

		static Scalar get(const PetscVector &v, const PetscInt index);
		static Scalar get(const PetscMatrix &v, const PetscInt row, const PetscInt col);
		static void get(const PetscVector &v, const std::vector<PetscInt> &index, std::vector<PetscScalar> &values);

		template<typename I>
		inline static void get(const PetscVector &v, const std::vector<I> &index, std::vector<PetscScalar> &values)
		{
			get(v, convert_to_petsc(index), values);
		}

		template<typename I>
		inline static std::vector<PetscInt> convert_to_petsc(const std::vector<I> &index)
		{
			std::vector<PetscInt> petsc_index(index.size());
			std::copy(index.begin(), index.end(), petsc_index.begin());
			return petsc_index;
		}

		inline static const std::vector<PetscInt> & convert_to_petsc(const std::vector<PetscInt> &index)
		{
			return index;
		}

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
			result = v;
			result.scale(-1.);
		}

		//[binary]
		static void apply_binary(PetscVector &result, const Reciprocal<Scalar> &reciprocal, const PetscVector &vec);
		static void apply_binary(Vector &result, const Scalar factor, const Multiplies &, const Vector &vec);
		static void apply_binary(Matrix &result, const Scalar factor, const Multiplies &, const Matrix &mat);
		static void apply_binary(PetscVector &result, const PetscMatrix &left, const Multiplies &, const PetscVector &right);
		static void apply_binary(PetscMatrix &result, const PetscMatrix &left, const Multiplies &, const PetscMatrix &right);
		static void apply_binary(PetscVector &result, const PetscVector &left, const EMultiplies &, const PetscVector &right);
		static void apply_binary(PetscVector &result, const PetscVector &left, const Divides &, const PetscVector &right);
		static void apply_binary(PetscVector &result, const PetscVector &left, const Min &, const PetscVector &right);
		static void apply_binary(PetscVector &result, const PetscVector &left, const Max &, const PetscVector &right);

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

		static void allocate_apply_vec(PetscVector &result, const PetscVector &left, const PetscVector &right);

		template<class Operation>
		static void apply_binary_generic(PetscVector &result, const PetscVector &left, const Operation &op, const PetscVector &right)
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
		static void mat_mult_add(PetscVector &result, const PetscMatrix &m, const PetscVector &right, const PetscVector &left);
		static void mat_mult_t_add(PetscVector &result, const PetscMatrix &m, const PetscVector &right, const PetscVector &left);		
		static void triple_product_ptap(PetscMatrix &result, const PetscMatrix &, const PetscMatrix &); 
		static void triple_product(PetscMatrix &result, const PetscMatrix &, const PetscMatrix &, const PetscMatrix &); 


		static Scalar norm2(const PetscVector &v);
		static Scalar norm2(const Matrix &m);
		static Scalar norm1(const PetscVector &v);
		static Scalar norm_infty(const PetscVector &v);
		static Scalar reduce(const PetscVector &vec, const Plus &);
		static Scalar reduce(const PetscMatrix &mat, const Plus &op);
		static Scalar reduce(const PetscVector &, const Min &);
		static Scalar reduce(const PetscMatrix &, const Min &);
		static Scalar reduce(const PetscVector &, const Max &);
		static Scalar reduce(const PetscMatrix &, const Max &);
		static Scalar dot(const Vector &left, const Vector &right);
		
		// get diagonal of matrix as vector
		static void diag(PetscVector &vec, const PetscMatrix &mat);
		static void diag(PetscSparseMatrix &mat, const PetscVector &vec);
		static void diag(PetscMatrix &mat, const PetscVector &vec);
		static void diag(PetscMatrix &out, const PetscMatrix &in);
        static void mat_diag_shift(PetscMatrix &left, const Scalar diag_factor);
        

		static bool compare(const Vector &left, const Vector &right, const ApproxEqual &comp);
		static bool compare(const Matrix &left, const Matrix &right, const ApproxEqual &comp);
		
		
		static void kronecker_product(Matrix &result, const Vector &left, const Vector &right);
		
		//[selection]
		static void assign_from_range(PetscVector &left, const PetscVector &right, const Range &globalRowRange,
							 const Range &global_col_range);
		
		static void assign_to_range(PetscMatrix &left, const PetscMatrix &right, const Range &global_row_range,
						   const Range &global_col_range);
		
		static void assign_to_range( PetscMatrix &left, const Identity &, const Range &global_row_range,
						   const Range &global_col_range);
		
		static void assign_to_range( PetscVector &left, const PetscVector &right, const Range &global_row_range,
						   const Range &global_col_range);

		static void assign_from_range(
			PetscMatrix &left,
			const PetscMatrix &right,
			const Range &globalRowRange,
			const Range &global_col_range);

	

		static void select(
					PetscVector &left,
					const PetscVector &right,
		      		const std::vector<PetscInt> &index);

		static void select(
					PetscMatrix &left,
					const PetscMatrix &right,
		      		const std::vector<PetscInt> &row_index,
		      		const std::vector<PetscInt> &col_index);


	
		// this is not very correct yet ...
		static void build_local_diag_block(PetscSerialSparseMatrix &left, const PetscSparseMatrix &right);
		static void build_local_redistribute(PetscVector &result, const PetscVector &, const PetscVector &); 

		static bool is_nan_or_inf(const PetscVector &v); 
		static bool is_nan_or_inf(const PetscMatrix &m); 


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

		static void axpy(PetscVector &y, const Scalar alpha, const PetscVector &x);
		static void axpy(PetscMatrix &y, const Scalar alpha, const PetscMatrix &x);
		static void axpby(Vector &y, const Scalar alpha, const Vector &x, const Scalar &beta);

		inline static void multiply(
			PetscMatrix &result,
			bool transpose_left,
			const PetscMatrix &left,
			bool transpose_right,
			const PetscMatrix &right)
		{
			gemm(result, 0.0, 1., transpose_left, left, transpose_right, right);
		}

		inline static void multiply(
			Vector &result,
			bool transpose_left,
			const PetscMatrix &left,
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
			return PetscErrorHandler::Check(err);
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

		static void build_ghosts(
			const PetscInt &local_size,
			const PetscInt &size,
			const std::vector<PetscInt> &index,
			PetscVector &vec);

		template<typename I>
		inline static void build_ghosts(
			const PetscInt &local_size,
			const PetscInt &size,
			const std::vector<I> &index,
			PetscVector &vec)
		{
			build_ghosts(local_size, size, convert_to_petsc(index), vec);
		}


		static void update_ghosts(PetscVector &vec);

	private:

		static void gemm(
			PetscMatrix &result,
			const Scalar beta,
			const Scalar alpha,
			bool transpose_left,
			const PetscMatrix &left,
			bool transpose_right,
			const PetscMatrix &right);

		static void gemv(
			Vector &y, 
			const Scalar beta,
			const Scalar alpha,
			bool transpose_A,
			const PetscMatrix &A, 
			const Vector &x);
				

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

