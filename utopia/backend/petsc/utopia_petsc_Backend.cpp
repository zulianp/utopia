
#include "utopia_petsc_Backend.hpp"

#include "petscmat.h"
#include "petscvec.h"

typedef utopia::PetscBackend::Scalar Scalar;

namespace utopia {

	class CompatibleMatPair {
	public:
		CompatibleMatPair(const MPI_Comm comm, const Mat &left, const Mat &right)
		{
			const bool left_is_sparse  = is_sparse(left);
			const bool right_is_sparse = is_sparse(right);

			must_destroy_left_  = false;
			must_destroy_right_ = false;

			MatType common_type;
			if(left_is_sparse != right_is_sparse) {
				if(left_is_sparse) {
					MatCreate(comm, &left_);
					PetscBackend::check_error( MatGetType(right, &common_type) );
					PetscBackend::check_error( MatConvert(left, common_type, MAT_INITIAL_MATRIX, &left_) );
					right_ = right;
					must_destroy_left_  = true;

				} else {
					MatCreate(comm, &right_);
					PetscBackend::check_error( MatGetType(left, &common_type) );
					PetscBackend::check_error( MatConvert(right, common_type, MAT_INITIAL_MATRIX, &right_) );
					left_ = left;
					must_destroy_right_ = true;
				}
			} else {
				left_  = left;
				right_ = right;
			}
		}


		inline const Mat &left()
		{
			return left_;
		}

		inline const Mat &right()
		{
			return right_;
		}

		~CompatibleMatPair()
		{
			if(must_destroy_left_) MatDestroy(&left_);
			if(must_destroy_right_) MatDestroy(&right_);
		}

	private:

		static bool is_sparse(const Mat &mat)
		{
			MatType type;
			MatGetType(mat, &type);
			const std::string type_str(type);
			const size_t start = type_str.size() - 5;
			return !(type_str.substr(start, 5) == "dense");
		}

		Mat left_;
		Mat right_;
		bool must_destroy_left_;
		bool must_destroy_right_;
	};

	void PetscBackend::resize(PetscMatrix &mat, const Size &s)
	{
		Mat &m = mat.implementation();
		PetscInt r, c;
		MatGetSize(mat.implementation(), &r, &c);
		if(r == s.get(0) && c == s.get(1)) {
			return;
		}

		MPI_Comm comm = PetscObjectComm((PetscObject) mat.implementation());
		MatType type;
		check_error( MatGetType(m, &type) );
		MatDestroy(&m);

		MatCreate(comm, &m);
		MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1));
		MatSetType(m, type);
	}

	void PetscBackend::resize(PetscMatrix &mat, const Size &local_s, const Size &global_s)
	{
		Mat &m = mat.implementation();
		PetscInt r, c, R, C;
		MatGetLocalSize(mat.implementation(), &r, &c);
		MatGetSize(mat.implementation(), &R, &C);
		if(r == local_s.get(0) && c == local_s.get(1) && R == global_s.get(0) && C == global_s.get(1)) {
			return;
		}

		MPI_Comm comm = PetscObjectComm((PetscObject) mat.implementation());
		MatType type;
		check_error( MatGetType(m, &type) );
		MatDestroy(&m);

		MatCreate(comm, &m);
		MatSetSizes(m, local_s.get(0), local_s.get(1), global_s.get(0), global_s.get(1));
		MatSetType(m, type);
	}

	void PetscBackend::resize(PetscVector &vec, const Size &s)
	{
		Vec &v = vec.implementation();
		PetscInt n;

		if(vec.initialized()) {
			VecGetSize(v,  &n);
			
			if(n == s.get(0)) {
				return;
			}
		}

		MPI_Comm comm; 
		if(!vec.is_null()) {
			comm = PetscObjectComm((PetscObject) vec.implementation());
			VecDestroy(&v);	
		} else {
			comm = default_communicator();
		}

		vec_create_parallel(comm, PETSC_DECIDE, s.get(0), &v);
	}

	void PetscBackend::resize(PetscVector &vec, const Size &s_local, const Size &s_global)
	{
		Vec &v = vec.implementation();
		PetscInt n, N;

		if(vec.initialized()) {
			VecSetSizes(v, s_local.get(0), s_global.get(0));
			return;
		}

		MPI_Comm comm; 
		if(!vec.is_null()) {
			comm = PetscObjectComm((PetscObject) vec.implementation());
			VecDestroy(&v);	
		} else {
			comm = default_communicator();
		}
		
		vec_create_parallel(comm, s_local.get(0), s_global.get(0), &v);
	}

	void PetscBackend::assign_transposed(PetscMatrix &left, const PetscMatrix &right)
	{
		right.transpose(left);
	}
	
	void PetscBackend::clear(PetscMatrix &mat)
	{
		mat.clear();
	}
	
	void PetscBackend::convert(Vec vec, PetscVector &wrapper)
	{
		wrapper.copy_from(vec);
	}
	
	void PetscBackend::convert(Mat mat, PetscMatrix &wrapper)
	{
		wrapper.copy_from(mat);
	}
	
	void PetscBackend::convert(Mat mat, PetscSparseMatrix &wrapper)
	{
		MatDestroy(&wrapper.implementation());
		MatDuplicate(mat, MAT_COPY_VALUES, &wrapper.implementation());
	}
	
	void PetscBackend::convert(const PetscVector &wrapper, Vec vec)
	{
		VecCopy(wrapper.implementation(), vec);
	}
	
	void PetscBackend::convert(const PetscMatrix &wrapper, Mat mat)
	{
		MatCopy(wrapper.implementation(), mat, DIFFERENT_NONZERO_PATTERN);
	}
	
	void PetscBackend::convert(const PetscSparseMatrix &wrapper, Mat mat)
	{
		MatCopy(wrapper.implementation(), mat, DIFFERENT_NONZERO_PATTERN);
	}

	void PetscBackend::wrap(Mat mat, PetscSparseMatrix &wrapper)
	{
		wrapper.wrap(mat);
	}

	void PetscBackend::wrap(Vec vec, PetscVector &wrapper)
	{
		assert(false && "TODO");
		// wrapper.wrap(vec);
	}

	Range PetscBackend::range(const PetscVector &v) 
	{
		return v.range();
	}
	
	Range PetscBackend::row_range(const PetscMatrix &m) 
	{
		return m.row_range();
	}
		

	//FIXME reinterpret meaning in dsl
	Range PetscBackend::col_range(const PetscMatrix &m) 
	{
		//FIXME do not use col_range to iterate
		// PetscInt grows, gcols;
		// MatGetSize(m.implementation(), &grows, &gcols);

		auto s = m.size();
		return Range(0, s.get(1));
	}
	
	void PetscBackend::size(const PetscMatrix &m, Size &size) 
	{
		size = m.size();
	}
	
	void PetscBackend::size(const PetscVector &v, Size &size) 
	{
		size.set_dims(1);
		size.set(0, v.size());
	}
	
	void PetscBackend::local_size(const PetscVector &v, Size &size) 
	{
		size.set_dims(1);
		size.set(0, v.local_size());
	}

	void PetscBackend::local_size(const PetscMatrix &mat, Size &size)
	{
		size = mat.local_size();
	}

	void PetscBackend::select(
		PetscVector &left,
		const PetscVector &right,
		const std::vector<PetscInt> &index)
	{
		right.select(index, left);
	}


	void PetscBackend::select(PetscMatrix &left,
		const PetscMatrix &right,
		const std::vector<PetscInt> &row_index,
		const std::vector<PetscInt> &col_index)
	{
		right.select(row_index, col_index, left);
	}

	void PetscBackend::assign_from_range(
		PetscMatrix &left,
		const PetscMatrix &right,
		const Range &global_row_range,
		const Range &global_col_range) 
	{
		assert(!global_row_range.empty());
		right.select(global_row_range, global_col_range, left);
	}
	
	// read matrix
	bool PetscBackend::read(const std::string &path, PetscMatrix &mat, const PetscArgs &args)
	{
		 if(!mat.read(args.comm, path)) {
		 	return false;
		 }

		 apply_args(args, mat);
		 return true;
	}

	bool is_matlab_file(const std::string &path)
	{
		size_t pos = path.find_last_of(".");

		bool is_matlab = false;
		if(pos != std::string::npos) {
			is_matlab = path.substr(pos, path.size()-pos) == ".m";
		}

		return is_matlab;
	}
	
	// write matrix
	bool PetscBackend::write(const std::string &path, const PetscMatrix &m)
	{
		if(is_matlab_file(path)) { 
			return m.write_matlab(path);
		} else { 
			return m.write(path);
		}
	}
	
	// write vector
	bool PetscBackend::write(const std::string &path, const PetscVector &v)
	{
		if(is_matlab_file(path)) { 
			return v.write_matlab(path);
		} else { 
			return v.write(path);
		}
	}
	
	void PetscBackend::monitor(const long &iteration, PetscMatrix &m)
	{
		PetscViewer viewer = nullptr;

		PetscViewerASCIIOpen(m.communicator(), "log_hessian.m", &viewer);  
		PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); 

		const char *name;
		PetscObjectGetName((PetscObject)m.implementation(), &name);
		PetscObjectSetName((PetscObject)m.implementation(), ("H_" + std::to_string(iteration)).c_str());
		MatView(m.implementation(), viewer); 

		PetscViewerDestroy(&viewer);
		PetscObjectSetName((PetscObject)m.implementation(), name);
	}
	
	void PetscBackend::monitor(const long &iteration, PetscVector &v)
	{
		PetscViewer viewer = nullptr;
        
		PetscViewerASCIIOpen(v.communicator(), "log_iterate.m", &viewer);  
		PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

		const char *name;
		PetscObjectGetName((PetscObject)v.implementation(), &name);
		PetscObjectSetName((PetscObject)v.implementation(), ("it_" + std::to_string(iteration)).c_str());
		VecView(v.implementation(), viewer); 
		PetscViewerDestroy(&viewer);
		PetscObjectSetName((PetscObject)v.implementation(), name);
	}

	Scalar PetscBackend::get_global_nnz(PetscMatrix &m)
	{	
		MatInfo        info;
		MatGetInfo(m.implementation(), MAT_GLOBAL_SUM, &info); 
		return info.nz_used; 
	}

	Scalar PetscBackend::get_local_nnz(PetscMatrix &m)
	{	
		MatInfo        info;
		MatGetInfo(m.implementation(), MAT_LOCAL, &info); 
		return info.nz_used; 
	}

	void PetscBackend::set_zero_rows(PetscMatrix &Mat_A, const std::vector<int> &index)
	{
		check_error(MatZeroRows(Mat_A.implementation(), index.size(), &index[0], 1.0, nullptr, nullptr));
	}

	void PetscBackend::apply_BC_to_system(PetscMatrix & A, PetscVector& x, PetscVector& rhs, const std::vector<int> &index)
	{
		check_error(MatZeroRows(A.implementation(), index.size(), &index[0], 1.0, x.implementation(), rhs.implementation())); 
	}

	// read vector
	bool PetscBackend::read(const std::string &path, PetscVector &vec, const PetscArgs &args)
	{
		if(!vec.read(args.comm, path)) return false;
		apply_args(args, vec);
		return true;
	}
	
	void PetscBackend::assign_from_range(
		PetscVector &left, 
		const PetscVector &right, 
		const Range &global_row_range,
		const Range &) 
	{
		right.select(global_row_range, left);
	}
	
	//remove me and make blas specialization
	void PetscBackend::gemm(
		PetscMatrix &result,
		const Scalar beta,
		const Scalar alpha,
		bool transpose_left,
		const PetscMatrix &left,
		bool transpose_right,
		const PetscMatrix &right) 
	{

		//TODO


		//FIXME only works for beta == 0 for the moment
		assert(fabs(beta) < 1e-16);
		
		//FIXME only works for alpha == 1 for the moment
		assert(fabs(alpha - 1) < 1e-16);

		CompatibleMatPair mat_pair(left.communicator(), left.implementation(), right.implementation());
		auto l = mat_pair.left();
		auto r = mat_pair.right();

		// MatDestroy( &result.implementation());
		result.destroy();

		bool ok = false;
		if(transpose_left && !transpose_right) {
			ok = check_error( MatTransposeMatMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
		} else if(!transpose_left && transpose_right) {
			ok = check_error( MatMatTransposeMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
		} else if(!transpose_left && !transpose_right) {
			ok = check_error( MatMatMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
		} else {
			assert(transpose_left && transpose_right);
			PetscMatrix temp;

			if(!check_error( MatMatMult(r, l, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation()) ) ) {
				ok = false;
			} else {
				assign_transposed(result, temp);
				ok = true;
			}
		}

		assert(ok);
	}

	//remove me and make blas specialization
	void PetscBackend::gemv(
		Vector &y, 
		const Scalar beta,
		const Scalar alpha,
		bool transpose_A,
		const PetscMatrix &A, 
		const Vector &x)
	{
		assert(fabs(beta) < 1e-16);
		assert(fabs(alpha - 1.) < 1e-16);
		
		if(transpose_A) {
			A.mult_t(x, y);
		} else {
			A.mult(x, y);
		}
		
	}

	void PetscBackend::build(PetscMatrix &m, const Size &size, const Identity &, const PetscArgs &opts)
	{
		m.dense_init_identity(opts.comm, parallel_dense_matrix_type(), PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), 1.);
		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const Identity &, const PetscArgs &opts)
	{
		m.matij_init_identity(opts.comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), 1.);
		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts)
	{
		m.dense_init_identity(opts.comm, parallel_dense_matrix_type(), size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, 1.);
		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts)
	{
		m.matij_init_identity(opts.comm, size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, 1.);
		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz, const PetscArgs &opts)
	{
		m.matij_init(
        	opts.comm,
        	PETSC_DECIDE,
        	PETSC_DECIDE,
        	size.get(0),
        	size.get(1),
        	nnz.nnz(),
        	nnz.nnz()
        );

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz, const PetscArgs &opts)
	{
		m.matij_init(
        	opts.comm,
        	size.get(0),
        	size.get(1),
        	PETSC_DETERMINE,
        	PETSC_DETERMINE,
        	nnz.nnz(),
        	nnz.nnz()
        );

		apply_args(opts, m);
	}
	
	/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
	void PetscBackend::build(PetscMatrix  &m, const Size &size, const LocalNNZ<PetscInt> &, const PetscArgs &opts)
	{
		build(m, size, LocalValues<Scalar>(0));
	}
	
	/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
	void PetscBackend::build(PetscMatrix  &m, const Size &size, const NNZ<PetscInt> &, const PetscArgs &opts)
	{
		build(m, size, Zeros());
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const Zeros &, const PetscArgs &opts)
	{
		build(m, size, Values<Scalar>(0));
	}
	
	void PetscBackend::build(PetscVector &v, const Size &size, const Zeros &, const PetscArgs &opts)
	{
		build(v, size, Values<Scalar>(0));
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const LocalZeros &, const PetscArgs &opts)
	{
		build(m, size, LocalValues<Scalar>(0));
	}
	
	void PetscBackend::build(PetscVector &v, const Size &size, const LocalZeros &, const PetscArgs &opts)
	{
		build(v, size, LocalValues<Scalar>(0));
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const Values<Scalar> &values, const PetscArgs &opts) {
		m.dense_init_values(opts.comm, parallel_dense_matrix_type(), PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), values.value());
		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts) {
		m.dense_init_values(opts.comm, parallel_dense_matrix_type(), size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, values.value());
		apply_args(opts, m);
	}

	void PetscBackend::build(PetscVector &v, const Size &in_local_size, const Size &&in_global_size, const Values<Scalar> &values, const PetscArgs &opts)
	{
		v.values(opts.comm, parallel_vector_type(), in_local_size.get(0), in_global_size.get(0), values.value());
		apply_args(opts, v);
	}
	
	void PetscBackend::build(PetscVector &v, const Size &size, const Values<Scalar> &values, const PetscArgs &opts)
	{
		v.values(opts.comm, parallel_vector_type(), PETSC_DECIDE, size.get(0), values.value());
		apply_args(opts, v);
	}
	
	void PetscBackend::build(PetscVector &v, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts)
	{
		v.values(opts.comm, parallel_vector_type(), size.get(0), PETSC_DETERMINE, values.value());
		apply_args(opts, v);
	}
	
	void PetscBackend::add(PetscVector &v, const PetscInt index, Scalar value)
	{
		v.add(index, value);
	}
	
	void PetscBackend::add(PetscMatrix &m, const PetscInt row, const PetscInt col, Scalar value)
	{
		m.add(row, col, value);
	}
	
	void PetscBackend::set(PetscVector &v, const PetscInt index, Scalar value)
	{
		v.set(index, value);
	}
	
	void PetscBackend::set(PetscVector &v, const std::vector<PetscInt> &indices, const std::vector<Scalar> &values)
	{	
		v.set(indices, values);
	}
	
	void PetscBackend::set(PetscMatrix &m, const PetscInt row, const PetscInt col, Scalar value)
	{
		m.set(row, col, value);
	}
	
	void PetscBackend::set(PetscSparseMatrix &m, const PetscInt row, const PetscInt col, Scalar value)
	{
		m.set(row, col, value);
	}
	
	void PetscBackend::write_lock(PetscVector &vec)
	{
		vec.write_lock();
	}
	
	void PetscBackend::write_unlock(PetscVector &vec)
	{
		vec.write_unlock();
	}
	
	void PetscBackend::write_lock(PetscMatrix &mat) 
	{
		mat.write_lock();
	}
	
	void PetscBackend::write_unlock(PetscMatrix &mat)
	{
		mat.write_unlock();
	}
	
	//FIXME remove method from dsl
	void PetscBackend::set(
		PetscMatrix &m,
		const std::vector<PetscInt> &rows,
		const std::vector<PetscInt> &cols,
		const std::vector<Scalar> &values) 
	{
		assert(rows.size() == values.size());
		assert(cols.size() == values.size());

		for(std::vector<PetscInt>::size_type i = 0; i != rows.size(); ++i) {
			set(m, rows[i], cols[i], values[i]);
		}
	}

	void PetscBackend::add_matrix(
		PetscMatrix &m,
		const std::vector<PetscInt> &rows,
		const std::vector<PetscInt> &cols,
		const std::vector<Scalar> &values)
	{
		m.add_matrix(rows, cols, values);
	}

	Scalar PetscBackend::get(const PetscVector &v, const PetscInt index) {
		return v.get(index);
	}
	
	Scalar PetscBackend::get(const PetscMatrix &m, const PetscInt row, const PetscInt col) {
		return m.get(row, col);
	}

	void PetscBackend::get(const PetscVector &v, const std::vector<PetscInt> &index, std::vector<PetscScalar> &values)
	{
		v.get(index, values);
	}
	
	void PetscBackend::apply_binary(PetscVector &result, const PetscMatrix &left, const Multiplies &, const PetscVector &right)
	{
		left.mult(right, result);
	}
	
	void PetscBackend::apply_binary(PetscMatrix &result, const PetscMatrix &left, const Multiplies &, const PetscMatrix &right)
	{
		left.mult(right, result);
	}

	void PetscBackend::mat_mult_add(PetscVector &result, const PetscMatrix &m, const PetscVector &v1, const PetscVector &v2)
	{
		m.mult_add(v1, v2, result);
	}

	void PetscBackend::mat_mult_t_add(PetscVector &result, const PetscMatrix &m, const PetscVector &v1, const PetscVector &v2)
	{
		m.mult_t_add(v1, v2, result);
	}

	void PetscBackend::apply_unary(Vector &result, const Abs &, const Vector &vec)
	{
		result = vec;
		result.abs();
	}

	void PetscBackend::allocate_apply_vec(PetscVector &result, const PetscVector &left, const PetscVector &right)
	{
		if(result.implementation() != left.implementation() && result.implementation() != right.implementation()) {
			vec_repurpose(
				right.communicator(),
				right.type(),
				right.local_size(),
				right.size(),
				&result.implementation()
			);
		}
	}
	
	void PetscBackend::apply_binary(PetscVector &result, const PetscVector &left, const EMultiplies &, const PetscVector &right)
	{
		left.e_mul(right, result);
	}
	
	void PetscBackend::apply_binary(PetscVector &result, const PetscVector &left, const Divides &, const PetscVector &right)
	{
		left.e_div(right, result);
	}

	void PetscBackend::apply_binary(PetscVector &result, const PetscVector &left, const Min &op, const PetscVector &right)
	{
		apply_binary_generic(result, left, op, right);
	}

	void PetscBackend::apply_binary(PetscVector &result, const PetscVector &left, const Max &op, const PetscVector &right)
	{
		apply_binary_generic(result, left, op, right);
	}

	// reciprocal
	void PetscBackend::apply_binary(PetscVector &result, const Reciprocal<Scalar> &reciprocal, const PetscVector &vec)
	{
		result = vec;
		result.reciprocal(reciprocal.numerator());
	}
	
	Scalar PetscBackend::norm2(const PetscVector &v)
	{
		return v.norm2();
	}
	
	Scalar PetscBackend::norm2(const Matrix &m)
	{
		return m.norm2();
	}
	
	Scalar PetscBackend::norm1(const PetscVector &v) {
		return v.norm1();
	}
	
	Scalar PetscBackend::norm_infty(const PetscVector &v) {
		return v.norm_infty();
	}
	
	Scalar PetscBackend::reduce(const PetscVector &vec, const Plus &) {
		return vec.sum();
	}
	
	Scalar PetscBackend::reduce(const PetscMatrix &mat, const Plus &op) {
		return mat.sum();
	}
	
	Scalar PetscBackend::reduce(const PetscVector &v, const Min &)
	{
		return v.min();
	}

	Scalar PetscBackend::reduce(const PetscMatrix &m, const Min &op)
	{
		return m.min();
	}

	Scalar PetscBackend::reduce(const PetscVector &v, const Max &)
	{
		return v.max();
	}

	Scalar PetscBackend::reduce(const PetscMatrix &m, const Max &op)
	{
		return m.max();
	}

	// get diagonal of matrix as vector
	void PetscBackend::diag(PetscVector &vec, const PetscMatrix &mat)
	{
		mat.get_diag(vec);
	}
	
	void PetscBackend::diag(PetscSparseMatrix &mat, const PetscVector &vec)
	{
		//FIXME make it generic for all matrix types
		mat.matij_init_diag(vec);
	}

	void PetscBackend::diag(PetscMatrix &out, const PetscMatrix &in)
	{
		in.get_diag(out);
	}

	void PetscBackend::diag(PetscMatrix &mat, const PetscVector &vec)
	{
		//FIXME make it generic for all matrix types
		mat.dense_init_diag(parallel_dense_matrix_type(), vec);
	}

	void PetscBackend::mat_diag_shift(PetscMatrix &left, const Scalar diag_factor)
	{
		check_error( MatShift(left.implementation(), diag_factor) );
	}
	
	bool PetscBackend::compare(const Vector &left, const Vector &right, const ApproxEqual &comp) {
		PetscVector diff;
		apply_binary(diff, left, Minus(), right);
		return norm_infty(diff) <= comp.tol();
	}
	
	bool PetscBackend::compare(const Matrix &left, const Matrix &right, const ApproxEqual &comp) {
		Matrix diff;
		apply_binary(diff, left, Minus(), right);
		return norm2(diff) <= comp.tol();
	}
	
	void PetscBackend::apply_binary(Vector &result, const Scalar factor, const Multiplies &, const Vector &vec) {
		result = vec;
		result.scale(factor);
	}
	
	void PetscBackend::apply_binary(Matrix &result, const Scalar factor, const Multiplies &, const Matrix &mat) {
		result = mat;
		MatScale(result.implementation(), factor);
	}

	void PetscBackend::scale(Vector &result, const Scalar scale_factor)
	{
		result.scale(scale_factor);
	}

	void PetscBackend::scale(Matrix &result, const Scalar scale_factor)
	{
		MatScale(result.implementation(), scale_factor);
	}

	Scalar PetscBackend::dot(const Vector &left, const Vector &right) {
		return left.dot(right);
	}
	
	void PetscBackend::kronecker_product(Matrix &result, const Vector &left, const Vector &right) {
		MPI_Comm comm;
		PetscObjectGetComm((PetscObject)right.implementation(),&comm);
		
		int n_procs = 0, rank = 0;
		MPI_Comm_size(comm, &n_procs);
		MPI_Comm_rank(comm, &rank);
		
		Size lsize, rsize;
		size(left, lsize);
		size(right, rsize);
		const Size result_size = {lsize.get(0), rsize.get(0)};
		
		const Scalar * right_array = nullptr;
		PetscErrorCode err =  VecGetArrayRead(right.implementation(), &right_array);
		
		const Scalar * left_array = nullptr;
		err =  VecGetArrayRead(left.implementation(), &left_array);
		
		const Range r_range = range(right);
		const Range l_range = range(left);
		const PetscInt n = rsize.get(0);
		
		std::vector<Scalar> recvbuf(n, 0);
		
		//not very efficient but good enough for the moment
		int is_evenly_distributed = n == r_range.extent() * n_procs;
		MPI_Allreduce(MPI_IN_PLACE, &is_evenly_distributed, 1, MPI_INT, MPI_MIN, comm);
		
		if(is_evenly_distributed) {
			MPI_Allgather(right_array, r_range.extent(), MPIU_SCALAR, &recvbuf[0], r_range.extent(), MPIU_SCALAR, comm);
		} else {
			const long n_values = r_range.extent();
			std::vector<long> n_values_x_proc(n_procs);
			std::vector<long> offsets(n_procs);
			n_values_x_proc[rank] = n_values;
			
			MPI_Allgather(&n_values, 1, MPI_LONG, &n_values_x_proc[0], 1, MPI_LONG, comm);
			
			offsets[0] = 0;
			for(int r = 1; r != n_procs; ++r) {
				offsets[r] = offsets[r-1] + n_values_x_proc[r-1];
			}
			
			std::copy(right_array, right_array+n_values, recvbuf.begin()+offsets[rank]);
			std::vector<MPI_Request> requests((n_procs-1)*2);
			
			SizeType req_index = 0;
			for(int r = 0; r < n_procs; ++r) {
				if(r == rank) continue;
				
				MPI_Isend(right_array,  r_range.extent(), MPIU_SCALAR, r, r, comm, &requests[req_index++]);
				MPI_Irecv(&recvbuf[offsets[r]], n_values_x_proc[r], MPIU_SCALAR, r, rank, comm, &requests[req_index++]);
			}
			
			MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
		}
		
		err = VecRestoreArrayRead(right.implementation(), &right_array);
		
		MatDestroy(&result.implementation());

		dense_mat_create_parallel(
			comm,
			l_range.extent(),
			PETSC_DECIDE,
			result_size.get(0),
			result_size.get(1),
			&result.implementation()
		);
		
		write_lock(result);
		
		for(SizeType i = l_range.begin(); i != l_range.end(); ++i) {
			const Scalar l_value = left_array[i-l_range.begin()];
			
			for(SizeType j = 0; j != n; ++j) {
				const Scalar r_value = recvbuf.at(j);
				
				MatSetValue(result.implementation(), i, j, l_value*r_value, INSERT_VALUES);
			}
		}
		
		write_unlock(result);
		
		VecRestoreArrayRead(left.implementation(), &left_array);
	}
	

	void PetscBackend::assign_to_range(PetscMatrix & /*left*/, const PetscMatrix &/*right*/, const Range &/*global_row_range*/,
									 const Range &/*global_col_range*/)
	{
		
		assert(false); //TODO
		
	}
	
	void PetscBackend::assign_to_range( PetscMatrix &/*left*/, const Identity &/**/, const Range &/*global_row_range*/,
									 const Range &/*global_col_range*/) {
		assert(false); //TODO
	}
	
	void PetscBackend::assign_to_range( PetscVector &/*left*/, const PetscVector &/*right*/, const Range &/*global_row_range*/,
									 const Range &/*global_col_range*/) {
		assert(false); //TODO
	}
	
	void PetscBackend::vec_to_mat(Matrix &m, const Vector &v, const bool transpose) {
		Matrix temp_transpose;
		Matrix &work = transpose ? temp_transpose : m;

		MatDestroy(&m.implementation());
		MatDestroy(&work.implementation());

		dense_mat_create_parallel(v.communicator(), v.local_size(), 1, v.size(), 1, &work.implementation());
		
		Range r = range(v);
		for (auto i = r.begin(); i < r.end(); ++i) {
			set(work, i, 0, get(v, i));
		}
		
		MatAssemblyBegin(work.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(work.implementation(), MAT_FINAL_ASSEMBLY);
		
		if (transpose) {
			check_error( MatTranspose(temp_transpose.implementation(), MAT_INITIAL_MATRIX, &m.implementation()) );
		}

	}
	
	// this is not very correct yet ...
	void PetscBackend::build_local_diag_block(PetscSerialSparseMatrix &left, const PetscSparseMatrix &right)
	{
		Mat M;
		check_error( MatGetDiagonalBlock(right.implementation(), &M) );
		convert(M, left);
	}
	
	void PetscBackend::triple_product_ptap(PetscMatrix & result, const PetscMatrix & A, const PetscMatrix &P)
	{
		if(result.implementation() != A.implementation() && result.implementation() != P.implementation()) {
			MatDestroy(&result.implementation());
		} else {
			std::cerr << "[Error] not handled case in triple_product_ptap" << std::endl;
		}

		check_error( MatPtAP(A.implementation(), P.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) ); 
	}
	
	void PetscBackend::triple_product(PetscMatrix & result, const PetscMatrix &A, const PetscMatrix &B, const PetscMatrix &C)
	{
		if(result.implementation() != A.implementation() && result.implementation() != B.implementation() && result.implementation() != C.implementation()) {
			result.destroy();
		}
		
		check_error( MatMatMatMult(A.implementation(), B.implementation(), C.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
	}

	bool PetscBackend::is_nan_or_inf(const PetscVector &X)
	{
		return X.is_nan_or_inf();
	}

	bool PetscBackend::is_nan_or_inf(const PetscMatrix &m)
	{
		return m.is_nan_or_inf();
	}

	// redistribution of local sizes of vector
	// TODO: can be done also for matrices 
	//       can be done also based on local sizes, no provided vector
	void PetscBackend::build_local_redistribute(PetscVector &result, const PetscVector &x_from, const PetscVector &shape_vec)
	{
		if(mpi_world_size() == 1) 
		{	
			result = x_from; 
		} 
		else
		{               
			PetscVector x_to = shape_vec; 
			IS is; 
			VecScatter newctx; 

			PetscInt ed, st; 
			VecGetOwnershipRange(x_to.implementation() ,&st, &ed);
			ISCreateStride(PETSC_COMM_SELF, (ed - st), st, 1, &is);
			VecScatterCreate(x_from.implementation(), is, x_to.implementation(), is, &newctx);
			VecScatterBegin(newctx, x_from.implementation(), x_to.implementation(), INSERT_VALUES, SCATTER_FORWARD);  
			VecScatterEnd(newctx, x_from.implementation(), x_to.implementation(), INSERT_VALUES, SCATTER_FORWARD);  


			result = x_to; 

			VecScatterDestroy(&newctx);
			ISDestroy(&is);
		}	
	}

	void PetscBackend::axpy(PetscVector &y, const Scalar alpha, const PetscVector &x) 
	{
		y.axpy(alpha, x);
	}

	void PetscBackend::axpy(PetscMatrix &y, const Scalar alpha, const PetscMatrix &x)
	{
		y.axpy(alpha, x);
	}

	void PetscBackend::axpby(Vector &y, const Scalar alpha, const Vector &x, const Scalar &beta)
	{
		y.axpby(alpha, x, beta);
	}
	
	void PetscBackend::diag_scale_right(Matrix &result, const Matrix &m, const Vector &diag)
	{
		result = m;
		check_error( MatDiagonalScale(result.implementation(), nullptr, diag.implementation()) );
	}

	void PetscBackend::diag_scale_left(Matrix &result, const Vector &diag, const Matrix &m)
	{
		result = m;
		check_error( MatDiagonalScale(result.implementation(), diag.implementation(), nullptr) );
	}

	void PetscBackend::diag_scale_left(Vector &result, const Vector &diag, const Vector &m)
	{
		diag.e_mul(m, result);
	}

	Scalar PetscBackend::trace(const Matrix &mat)
	{
		return mat.trace();
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Plus &, const int dim)
	{
		if(dim == 1) {
			mat.row_sum(result);
		} else {
			assert(false && "not available in pestsc");
			result.set_initialized(false);
		}
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Min &op, const int dim)
	{
		if (dim == 1) {
			mat.row_min(result);
		} else {
			assert(false && "not available in petsc");
			result.set_initialized(false);
		}
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Max &op, const int dim)
	{
		if (dim == 1) {
			mat.row_max(result);
		} else {
			assert(false && "not available in petsc");
			result.set_initialized(false);
		}
	}

	void PetscBackend::inverse(Matrix &result, const Matrix &mat)
	{
		mat.inverse(result);
	}

	void PetscBackend::vec_create_parallel(MPI_Comm comm, PetscInt n_local, PetscInt n_global, Vec *vec)
	{
		check_error( VecCreate(comm, vec) );
		check_error( VecSetType(*vec, parallel_vector_type()) );
		check_error( VecSetSizes(*vec, n_local, n_global) );
	}

	void PetscBackend::vec_repurpose(
		MPI_Comm comm,
		VecType type,
		PetscInt n_local,
		PetscInt n_global,
		Vec *vec)
	{
		if(vec == nullptr) {
			VecCreate(comm, vec);
		} else {
			if(comm != PetscObjectComm((PetscObject)*vec)) {
				check_error( VecDestroy(vec) );
				check_error( VecCreate(comm, vec) );
			} else {
				PetscInt old_n_global;
				VecGetSize(*vec, &old_n_global);

				if(old_n_global == n_global) {
					PetscInt old_n_local;
					VecGetLocalSize(*vec, &old_n_local);
					assert(old_n_local == n_local && "We do not handle the local consistency. Explicitly set sizes in the initialization.");
					return;
				}
			}
		}
		
		check_error( VecSetType(*vec, type) );
		check_error( VecSetSizes(*vec, n_local, n_global) );
	}

	void PetscBackend::sparse_mat_create_parallel(
		MPI_Comm comm,
		PetscInt rows_local,
		PetscInt cols_local,
		PetscInt rows_global,
		PetscInt cols_global,
		PetscInt d_nnz,
		PetscInt o_nnz,
		Mat *mat)
	{
		check_error( MatCreate(comm, mat) );
		check_error( MatSetSizes(*mat, rows_local, cols_local, rows_global, cols_global) );
		
		check_error( MatSetType(*mat, parallel_sparse_matrix_type()) );
		check_error( MatSeqAIJSetPreallocation(*mat, PetscMax(d_nnz, 1), PETSC_NULL) );
		check_error( MatMPIAIJSetPreallocation(*mat, PetscMax(d_nnz, 1), PETSC_NULL, PetscMax(o_nnz, 1), PETSC_NULL) ); 

		check_error( MatSetOption(*mat, MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE) );
		check_error( MatSetOption(*mat, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) );
		check_error( MatSetOption(*mat, MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) );

		check_error( MatZeroEntries(*mat) );
	}

	void PetscBackend::dense_mat_create_parallel(
		MPI_Comm comm,
		PetscInt rows_local,
		PetscInt cols_local,
		PetscInt rows_global,
		PetscInt cols_global,
		Mat *mat)
	{
		check_error( MatCreate(comm, mat) );
		check_error( MatSetType(*mat, parallel_dense_matrix_type()) );
		check_error( MatSetSizes(*mat, rows_local, cols_local, rows_global, cols_global) );
		check_error( MatSetUp(*mat) );

		check_error( MatSetOption(*mat, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) );
		check_error( MatSetOption(*mat, MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) );
	}

	void PetscBackend::apply_args(const PetscArgs &args, Matrix &m)
	{
		if(!args.name.empty()) {
			m.set_name(args.name);
		}
	}

	void PetscBackend::apply_args(const PetscArgs &args, Vector &v)
	{
		if(!args.name.empty()) {
			v.set_name(args.name);
		}
	}

	void PetscBackend::build_ghosts(
		const PetscInt &local_size,
		const PetscInt &global_size,
		const std::vector<PetscInt> &index,
		PetscVector &vec)
	{
		vec.ghosted(default_communicator(), local_size, global_size, index);
	}

	void PetscBackend::update_ghosts(PetscVector &vec)
	{
		vec.update_ghosts();
	}

}

