
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

	void PetscBackend::assign_transposed(PetscMatrix &left, const PetscMatrix &right) {
		if(&left != &right) {
			MatDestroy(&left.implementation());
		} else {
			assert(false);
		}

		MatTranspose(right.implementation(), MAT_INITIAL_MATRIX, &left.implementation());
	}
	
	void PetscBackend::clear(PetscMatrix &mat)
	{
		MPI_Comm comm = mat.communicator();
		MatDestroy(&mat.implementation());
		MatCreate(comm, &mat.implementation());
	}
	
	void PetscBackend::convert(Vec vec, PetscVector &wrapper)
	{
		MPI_Comm comm = PetscObjectComm((PetscObject) vec);

		PetscInt n, local_n;
		VecGetSize(vec, &n);
		VecGetLocalSize(vec, &local_n);
		VecDestroy(&wrapper.implementation());
		vec_create_parallel(comm, local_n, n, &wrapper.implementation());
		VecCopy(vec, wrapper.implementation());
		wrapper.set_initialized(true);
	}
	
	void PetscBackend::convert(Mat mat, PetscMatrix &wrapper)
	{
		MPI_Comm comm = PetscObjectComm((PetscObject) mat);

		PetscInt r, c, local_r, local_c;
		MatGetSize(mat, &r, &c);
		MatGetLocalSize(mat, &local_r, &local_c);
		
		MatDestroy(&wrapper.implementation());
		MatCreateDense(comm, local_r, local_c, r, c, NULL,
			&wrapper.implementation());
		MatCopy(mat, wrapper.implementation(), SAME_NONZERO_PATTERN);
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
		PetscInt rbegin, rend;
		VecGetOwnershipRange(v.implementation(), &rbegin, &rend);
		return Range(rbegin, rend);
	}
	
	Range PetscBackend::row_range(const PetscMatrix &m) 
	{
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		assert(Range(rbegin, rend).valid());
		return Range(rbegin, rend);
	}
	
	Range PetscBackend::col_range(const PetscMatrix &m) 
	{
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		return Range(0, gcols);
	}
	
	void PetscBackend::size(const PetscMatrix &m, Size &size) 
	{
		PetscInt grows, gcols;
		size.set_dims(2);
		MatGetSize(m.implementation(), &grows, &gcols);
		size.set(0, grows);
		size.set(1, gcols);
	}
	
	void PetscBackend::size(const PetscVector &m, Size &size) 
	{
		PetscInt n;
		size.set_dims(1);
		VecGetSize(m.implementation(), &n);
		size.set(0, n);
	}
	
	void PetscBackend::local_size(const PetscVector &m, Size &size) 
	{
		PetscInt n;
		size.set_dims(1);
		VecGetLocalSize(m.implementation(), &n);
		size.set(0, n);
	}

	void PetscBackend::local_size(const PetscMatrix &mat, Size &size)
	{
		PetscInt n, m;
		size.set_dims(2);
		MatGetLocalSize(mat.implementation(), &n, &m);
		size.set(0, n);
		size.set(1, m);
	}

	void PetscBackend::select(
		PetscVector &left,
		const PetscVector &right,
		const std::vector<PetscInt> &index)
	{
		Vec &l = left.implementation();
		Vec r  = right.implementation();

		MPI_Comm comm = right.communicator();

		if(left.is_null()) {
			VecCreate(comm, &l);
		} else if(left.communicator() != right.communicator()) {
			VecDestroy(&l);
			VecCreate(comm, &l);
		}

		IS is_in;
		ISCreateGeneral(comm, index.size(), &index[0], PETSC_USE_POINTER, &is_in);
		VecScatter scatter_context;

		VecType type;
		VecGetType(r, &type);
		VecSetType(l, type);
		VecSetSizes(l, index.size(), PETSC_DETERMINE);

		VecScatterCreate(r, is_in, l, nullptr, &scatter_context);
		VecScatterBegin(scatter_context, r, l, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter_context, r, l, INSERT_VALUES, SCATTER_FORWARD);
		ISDestroy(&is_in);
		VecScatterDestroy(&scatter_context);
	}


	void PetscBackend::select(PetscMatrix &left,
		const PetscMatrix &right,
		const std::vector<PetscInt> &row_index,
		const std::vector<PetscInt> &col_index)
	{
		if(col_index.empty()) {
			PetscInt n_rows, n_cols;
			MatGetSize(right.implementation(), &n_rows, &n_cols);
			std::vector<PetscInt> all_index(n_cols);
			for(PetscInt i = 0; i < n_cols; ++i) {
				all_index[i] = i;
			}

			select_aux(left, right, row_index, all_index);

		} else {
			select_aux(left, right, row_index, col_index);
		}

	}

	void PetscBackend::select_aux(PetscMatrix &left,
		const PetscMatrix &right,
		const std::vector<PetscInt> &row_index,
		const std::vector<PetscInt> &col_index)
	{
		Mat r = right.implementation();

		PetscInt min_col = col_index[0], max_col = col_index[1];
		for(std::size_t i = 0; i < col_index.size(); ++i) {
			min_col = std::min(min_col, col_index[i]);
			max_col = std::max(max_col, col_index[i]);
		}

		MPI_Comm comm = PetscObjectComm((PetscObject)r);

		int vals[2] = { min_col, -max_col };
		MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_INT, MPI_MIN, comm);
		min_col =  vals[0];
		max_col = -vals[1];

		PetscInt global_cols = max_col - min_col + 1;
		PetscInt local_cols = PETSC_DECIDE;
		PetscSplitOwnership(comm, &local_cols, &global_cols);

		unsigned long offsets_in = row_index.size();
		unsigned long offset_out = 0;

		MPI_Exscan(
			&offsets_in,
			&offset_out,
			1, 
			MPI_UNSIGNED_LONG ,
			MPI_SUM,
			comm);

		offset_out += min_col;
		par_assign_from_local_is(
			row_index,
			col_index, 
			min_col,
			Range(offset_out, offset_out + local_cols), right, left);
	}


	void PetscBackend::par_assign_from_local_range(
		const Range &local_row_range,
		const Range &local_col_range,
		const Range &global_col_range,
		const PetscMatrix &right,
		PetscMatrix &result)
	{
		PetscErrorCode ierr = 0;

		Mat &l = result.implementation();
		const Mat r = right.implementation();

		std::vector<PetscInt> remote_rows;
		remote_rows.reserve(local_row_range.extent());
		for(PetscInt l_row = local_row_range.begin(); l_row < local_row_range.end(); ++l_row) {	
			remote_rows.push_back(l_row);
		}

		std::vector<PetscInt> remote_cols;
		remote_cols.reserve(global_col_range.extent());
		for(PetscInt l_col = global_col_range.begin(); l_col < global_col_range.end(); ++l_col) {	
			remote_cols.push_back(l_col);
		}

		par_assign_from_local_is(remote_rows, remote_cols, global_col_range.begin(), local_col_range, right, result);
	}


	void PetscBackend::par_assign_from_local_is(
		const std::vector<PetscInt> &remote_rows,
		const std::vector<PetscInt> &remote_cols,
		const PetscInt global_col_offset,
		const Range &local_col_range,
		const PetscMatrix &right,
		PetscMatrix &result)
	{

		Mat &l = result.implementation();
		const Mat r = right.implementation();


		std::stringstream ss;

		for(auto r : remote_rows) {
			ss << r;
		}

		IS isrow;
		PetscErrorCode ierr = ISCreateGeneral(right.communicator(), remote_rows.size(), &remote_rows[0], PETSC_USE_POINTER, &isrow);

		IS iscol;
		ierr = ISCreateGeneral(right.communicator(), remote_cols.size(), &remote_cols[0], PETSC_USE_POINTER, &iscol);

		MPI_Comm comm = PetscObjectComm((PetscObject)r);
		int size;
		MPI_Comm_size(comm, &size);

		int rank;
		MPI_Comm_rank(comm, &rank);

		//TODO maybe make it with collective comms for finiding out if there are off proc entries
		bool has_off_proc_entries = size > 1;

		if(has_off_proc_entries) {
			Mat * l_ptr;
			ierr = MatGetSubMatrices(r, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, &l_ptr);

			MatType type;
			MatGetType(r, &type);
			MatSetType(l, type);
			MatSetSizes(l, remote_rows.size(), local_col_range.extent(), PETSC_DETERMINE, PETSC_DETERMINE);
			MatSetUp(l);

			PetscInt rbegin, rend, cbegin, cend;
			MatGetOwnershipRange(l, &rbegin, &rend);
			MatGetOwnershipRangeColumn(l, &cbegin, &cend);
			PetscInt n_rows = rend - rbegin;

			PetscInt n_values;
			const PetscInt * cols;
			const Scalar * values;

			

			for(PetscInt row = 0; row < n_rows; ++row) {

				MatGetRow(*l_ptr, row, &n_values, &cols, &values);

				for(PetscInt i = 0; i < n_values; ++i) {
					MatSetValue(l, rbegin + row, cols[i] - global_col_offset, values[i], INSERT_VALUES);
				}

				MatRestoreRow(*l_ptr, row, &n_values, &cols, &values);
			}

			MatAssemblyBegin(l, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(l, MAT_FINAL_ASSEMBLY);

			MatDestroy(l_ptr);

		} else {
			MatDestroy(&l);
			ierr = MatGetSubMatrix(r, isrow, iscol, MAT_INITIAL_MATRIX, &l);
		}

		ISDestroy(&isrow);
		ISDestroy(&iscol);
	}

	void PetscBackend::par_assign_from_global_range(
		const Range &global_row_range,
		const Range &global_col_range,
		const PetscMatrix &right,
		PetscMatrix &result)
	{

		const Mat r = right.implementation();
		PetscInt global_rows = global_row_range.extent();
		PetscInt local_rows  = PETSC_DECIDE;

		PetscInt global_cols = global_col_range.extent();
		PetscInt local_cols  = PETSC_DECIDE;

		MPI_Comm comm = PetscObjectComm((PetscObject)r);
		PetscSplitOwnership(comm, &local_rows, &global_rows);
		PetscSplitOwnership(comm, &local_cols, &global_cols);

		unsigned long offsets_in[2] = {
			(unsigned long)local_rows,
			(unsigned long)local_cols
		};

		unsigned long offset_out[2] = {
			(unsigned long)0,
			(unsigned long)0
		};
		
		MPI_Exscan(
			&offsets_in,
			&offset_out,
			2, 
			MPI_UNSIGNED_LONG ,
			MPI_SUM,
			comm);

		offset_out[0] += global_row_range.begin();
		offset_out[1] += global_col_range.begin();

		par_assign_from_local_range(
			Range(offset_out[0], offset_out[0] + local_rows),
			Range(offset_out[1], offset_out[1] + local_cols),
			global_col_range,
			right,
			result);
	}
	
	void PetscBackend::assign_from_range(
		PetscMatrix &left,
		const PetscMatrix &right,
		const Range &global_row_range,
		const Range &global_col_range) 
	{
		assert(!global_row_range.empty());
		par_assign_from_global_range(
			global_row_range,
			global_col_range,
			right,
			left);
	}
	

	// read matrix
	bool PetscBackend::read(const std::string &path, PetscMatrix &mat, const PetscArgs &args)
	{
		MatDestroy(&mat.implementation());
		PetscViewer fd;
		PetscViewerBinaryOpen(args.comm, path.c_str(), FILE_MODE_READ, &fd);
		MatCreate(args.comm, &mat.implementation());
		const bool status = check_error( MatLoad(mat.implementation(), fd) );
		PetscViewerDestroy(&fd);

		apply_args(args, mat);
		return status;
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
		const Mat &A = m.implementation();

		const bool is_matlab = is_matlab_file(path);

		if(is_matlab) {
			PetscErrorCode ierr;
			PetscViewer fd;
			ierr = PetscViewerASCIIOpen(m.communicator(),path.c_str(), &fd); //CHKERRV(ierr);
			ierr = PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB); //CHKERRV(ierr);
			ierr = MatView(A, fd); //CHKERRV(ierr);
			PetscViewerDestroy(&fd);
			return check_error(ierr);
		} else {
			PetscViewer fd;
			PetscViewerBinaryOpen(m.communicator(), path.c_str(), FILE_MODE_WRITE, &fd);
			bool status;
			status =  check_error( MatView(A,fd));
			PetscViewerDestroy(&fd);
			return status;
		}
	}
	
	// write vector
	bool PetscBackend::write(const std::string &path, const PetscVector &v)
	{
		const Vec &A = v.implementation();

		bool is_matlab = is_matlab_file(path);
		if(is_matlab) {
			PetscErrorCode ierr;
			PetscViewer fd;
			ierr = PetscViewerASCIIOpen(v.communicator(), path.c_str(), &fd); //CHKERRV(ierr);
			ierr = PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB); //CHKERRV(ierr);
			ierr = VecView(A, fd); //CHKERRV(ierr);
			PetscViewerDestroy(&fd);
			return check_error(ierr);
		} else {
			PetscViewer fd;
			PetscViewerBinaryOpen(v.communicator(), path.c_str(), FILE_MODE_WRITE, &fd);
			bool status;
			status =  check_error( VecView(A,fd));
			PetscViewerDestroy(&fd);
			return status;
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
		VecDestroy(&vec.implementation());
		PetscViewer fd;
		PetscViewerBinaryOpen(args.comm, path.c_str(), FILE_MODE_READ, &fd);
		VecCreate(args.comm, &vec.implementation());
		const bool status = check_error( VecLoad(vec.implementation(), fd));
		PetscViewerDestroy(&fd);

		apply_args(args, vec);
		return status;
	}
	
	void PetscBackend::assign_from_range(
		PetscVector &left, 
		const PetscVector &right, 
		const Range &global_row_range,
		const Range &) 
	{
		assert(!global_row_range.empty());
		
		Range rr = range(right).intersect(global_row_range);
		
		vec_repurpose(
			right.communicator(),
			right.type(),
			PETSC_DECIDE,
			rr.extent(),
			&left.implementation()
		);

		write_lock(left);
		
		PetscInt r = 0;
		for (PetscInt rIt = rr.begin(); rIt < rr.end(); ++rIt) {
			r = rIt - global_row_range.begin();
			set(left, r, get(right, rIt));
		}
		
		write_unlock(left);
	}
	
	void PetscBackend::gemm(
		PetscMatrix &result,
		const Scalar beta,
		const Scalar alpha,
		bool transpose_left,
		const PetscMatrix &left,
		bool transpose_right,
		const PetscMatrix &right) 
	{
		//FIXME only works for beta == 0 for the moment
		assert(fabs(beta) < 1e-16);
		
		//FIXME only works for alpha == 1 for the moment
		assert(fabs(alpha - 1) < 1e-16);

		CompatibleMatPair mat_pair(left.communicator(), left.implementation(), right.implementation());
		auto l = mat_pair.left();
		auto r = mat_pair.right();

		MatDestroy( &result.implementation());

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

		const Mat &A_im = A.implementation();
		Vec &y_im = y.implementation();
		const Vec &x_im = x.implementation();

		if(y.is_null()) {
			VecCreate(x.communicator(), &y_im);
		} else if(y.communicator() != x.communicator()) {
			VecDestroy(&y_im);
			VecCreate(x.communicator(), &y_im);
		}
		
		Size gs, ls;
		size(A, gs);
		local_size(A, ls);

		VecType type;
		VecGetType(x_im, &type);
		VecSetType(y_im, type);

		if(transpose_A) {
			VecSetSizes(y_im, ls.get(1), gs.get(1));
			check_error( MatMultTranspose(A_im, x_im, y_im) );
		} else {
			VecSetSizes(y_im, ls.get(0), gs.get(0));
			check_error( MatMult(A_im, x_im, y_im) );
		}

		y.set_initialized(true);
	}

	void PetscBackend::build(PetscMatrix &m, const Size &size, const Identity &, const PetscArgs &opts) {
		MatDestroy(&m.implementation());

		//FIXME use this: MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],Scalar diag,Vec x,Vec b)
		MatCreateDense(opts.comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), NULL,
			&m.implementation());
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		rend = PetscMin(rend, gcols);
		for (PetscInt i = rbegin; i < rend; ++i) {
			MatSetValue(m.implementation(), i, i, 1, INSERT_VALUES);
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const Identity &, const PetscArgs &opts) {
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MatDestroy(&m.implementation());
		
		sparse_mat_create_parallel(
			opts.comm,
			PETSC_DECIDE,
			PETSC_DECIDE,
			rows,
			cols,
			1,
			1,
			&m.implementation());

		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		rend = PetscMin(rend, gcols);
		
		for (PetscInt i = rbegin; i < rend; ++i) {
			const PetscInt offset = i;
			const Scalar val = 1;
			MatSetValues(m.implementation(), 1, &offset, 1, &offset, &val, INSERT_VALUES);
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts) {
		MatDestroy(&m.implementation());

		//FIXME use this: MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],Scalar diag,Vec x,Vec b)
		MatCreateDense(opts.comm, size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, NULL,
			&m.implementation());
		
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows = size.get(0);
		PetscInt gcols = size.get(1);
		
		unsigned long send_buffer = gcols;
		unsigned long receive_buffer = 0;
		
		MPI_Exscan(&send_buffer, &receive_buffer, 1, MPI_UNSIGNED_LONG ,
			MPI_SUM, m.communicator());
		
		PetscInt extent = PetscMin(rend-rbegin, PetscMin(grows, gcols));
		
		for (PetscInt i = 0; i < extent; ++i) {
			MatSetValue(m.implementation(), rbegin+i, receive_buffer + i, 1, INSERT_VALUES);
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const LocalIdentity &, const PetscArgs &opts) {
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MatDestroy(&m.implementation());
		
		sparse_mat_create_parallel(
			opts.comm,
			rows,
			cols,
			PETSC_DETERMINE,
			PETSC_DETERMINE,
			1,
			1,
			&m.implementation());

		unsigned long send_buffer = cols;
		unsigned long receive_buffer = 0;
		
		MPI_Exscan(&send_buffer, &receive_buffer, 1, MPI_UNSIGNED_LONG ,
			MPI_SUM, m.communicator());
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt extent = PetscMin(rend-rbegin, PetscMin(rows, cols));

		check_error( MatZeroEntries(m.implementation()) );
		
		for (PetscInt i = 0; i < extent; ++i) {
			const PetscInt row_offset = rbegin+i;
			const PetscInt col_offset = receive_buffer+i;
			const Scalar val = 1;
			MatSetValues(m.implementation(), 1, &row_offset, 1, &col_offset, &val, INSERT_VALUES);
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz, const PetscArgs &opts) {
		
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);
		MatDestroy(&m.implementation());

		sparse_mat_create_parallel(
			opts.comm,
			PETSC_DECIDE,
			PETSC_DECIDE,
			rows,
			cols,
			nnz.nnz(),
			nnz.nnz(),
			&m.implementation());

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz, const PetscArgs &opts) {
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		sparse_mat_create_parallel(
			opts.comm,
			rows,
			cols,
			PETSC_DETERMINE,
			PETSC_DETERMINE,
			nnz.nnz(),
			nnz.nnz(),
			&m.implementation());

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
		MatDestroy(&m.implementation());

		MatCreateDense(opts.comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), NULL,
			&m.implementation());
		
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		check_error( MatZeroEntries(m.implementation()) );
		
		
		
		const Scalar v = values.value();
		for (PetscInt i = rbegin; i < rend; ++i) {
			for (PetscInt j = 0; j < gcols; ++j) {
				MatSetValue(m.implementation(), i, j, v, INSERT_VALUES);
			}
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscVector &v, const Size &in_local_size, const Size &&in_global_size, const Values<Scalar> &values, const PetscArgs &opts)
	{
		Size v_local_size;
		if(v.initialized()) {
			local_size(v, v_local_size);
		}	

		if(!v.initialized() || v_local_size.get(0) != in_local_size.get(0)) {
			VecDestroy(&v.implementation());
			vec_create_parallel(opts.comm, in_local_size.get(0), in_global_size.get(0), &v.implementation());
		}

		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());

		v.set_initialized(true);

		apply_args(opts, v);
	}
	
	void PetscBackend::build(PetscVector &v, const Size &size, const Values<Scalar> &values, const PetscArgs &opts) {

		VecDestroy(&v.implementation());
		vec_create_parallel(opts.comm, PETSC_DECIDE, size.get(0), &v.implementation());
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());

		v.set_initialized(true);

		apply_args(opts, v);
	}
	
	void PetscBackend::build(PetscMatrix &m, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		MatCreateDense(comm, size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, NULL,
			&m.implementation());
		
		//TODO check if it can be simplified using known information.
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		const Scalar v = values.value();
		for (PetscInt i = rbegin; i < rend; ++i) {
			for (PetscInt j = 0; j < gcols; ++j) {
				MatSetValue(m.implementation(), i, j, v, INSERT_VALUES);
			}
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);

		apply_args(opts, m);
	}
	
	void PetscBackend::build(PetscVector &v, const Size &size, const LocalValues<Scalar> &values, const PetscArgs &opts) {
		VecDestroy(&v.implementation());

		vec_create_parallel(opts.comm, size.get(0), PETSC_DETERMINE, &v.implementation());
		
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());

		v.set_initialized(true);

		apply_args(opts, v);
	}
	
	void PetscBackend::add(PetscVector &v, const PetscInt index, Scalar value)
	{
		VecSetValues(v.implementation(), 1, &index, &value, ADD_VALUES);
	}
	
	void PetscBackend::add(PetscMatrix &m, const PetscInt row, const PetscInt col, Scalar value)
	{
		MatSetValues(m.implementation(), 1, &row, 1, &col, &value, ADD_VALUES);
	}
	
	void PetscBackend::set(PetscVector &v, const PetscInt index, Scalar value) {
		VecSetValues(v.implementation(), 1, &index, &value, INSERT_VALUES);
	}
	
	void PetscBackend::set(PetscVector &v, const std::vector<PetscInt> &indices, const std::vector<Scalar> &values) {
		assert(indices.size() == values.size());
		VecSetValues(v.implementation(), indices.size(), &indices[0], &values[0], INSERT_VALUES);
	}
	
	
	void PetscBackend::set(PetscMatrix &v, const PetscInt row, const PetscInt col, Scalar value) {
		MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
	}
	
	void PetscBackend::set(PetscSparseMatrix &v, const PetscInt row, const PetscInt col, Scalar value) {
		MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
	}
	
	void PetscBackend::write_lock(PetscVector &) {}
	
	void PetscBackend::write_unlock(PetscVector &vec) {
		VecAssemblyBegin(vec.implementation());
		VecAssemblyEnd(vec.implementation());
		vec.set_initialized(true);
		vec.update_ghosts();
	}
	
	void PetscBackend::write_lock(const PetscMatrix &) { }
	
	void PetscBackend::write_unlock(const PetscMatrix &mat) {
		MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PetscBackend::write_lock(const PetscSparseMatrix &) { }
	
	void PetscBackend::write_unlock(const PetscSparseMatrix &mat) {
		MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
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
		assert(rows.size()*cols.size() == values.size());

		MatSetValues(
			m.implementation(), 
			static_cast<PetscInt>(rows.size()), &rows[0], 
			static_cast<PetscInt>(cols.size()), &cols[0],
			&values[0],
			ADD_VALUES);
	}

	Scalar PetscBackend::get(const PetscVector &v, const PetscInt index) {
		Scalar value;
		VecGetValues(v.implementation(), 1, &index, &value);
		return value;
	}
	
	Scalar PetscBackend::get(const PetscMatrix &v, const PetscInt row, const PetscInt col) {
		Scalar value;
		MatGetValues(v.implementation(), 1, &row, 1, &col, &value);
		return value;
	}

	void PetscBackend::get(const PetscVector &v, const std::vector<PetscInt> &index, std::vector<PetscScalar> &values)
	{
		v.get(index, values);
	}
	
	void PetscBackend::apply_binary(PetscVector &result, const PetscMatrix &left, const Multiplies &, const PetscVector &right) {
		auto sl = left.local_size();
		auto sg = left.size();

		vec_repurpose(right.communicator(), right.type(), sl.get(0), sg.get(0), &result.implementation());
		MatMult(left.implementation(), right.implementation(), result.implementation());
		result.set_initialized(true);
	}
	
	void PetscBackend::apply_binary(PetscMatrix &result, const PetscMatrix &left, const Multiplies &, const PetscMatrix &right)
	{
		if(right.implementation() != result.implementation() || left.implementation() != result.implementation()) {
			MatDestroy(&result.implementation());
			MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
		} else {
			PetscMatrix temp;
			MatDestroy(&temp.implementation());
			MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation());
			result = std::move(temp);
		}
	}

	void PetscBackend::mat_mult_add(PetscVector &result, const PetscMatrix &m, const PetscVector &v1, const PetscVector &v2)
	{
		if (v1.implementation() == result.implementation() || v2.implementation() == result.implementation()) {
			PetscVector temp;	
			vec_repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size(), &temp.implementation());
			MatMultAdd(m.implementation(), v1.implementation(), v2.implementation(), temp.implementation());
			result = std::move(temp);
		} else {
			vec_repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size(), &result.implementation());
			MatMultAdd(m.implementation(), v1.implementation(), v2.implementation(), result.implementation());
		}
	}

	void PetscBackend::mat_mult_t_add(PetscVector &result, const PetscMatrix &m, const PetscVector &v1, const PetscVector &v2)
	{
		if (v1.implementation() == result.implementation() || v2.implementation() == result.implementation()) {
			PetscVector temp;	
			vec_repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size(), &temp.implementation());
			MatMultTransposeAdd(m.implementation(), v1.implementation(), v2.implementation(), temp.implementation());
			result = std::move(temp);
		} else {
			vec_repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size(), &result.implementation());
			MatMultTransposeAdd(m.implementation(), v1.implementation(), v2.implementation(), result.implementation());
		}

		result.set_initialized(true);
	}

	void PetscBackend::apply_unary(Vector &result, const Abs &, const Vector &vec)
	{
		if(vec.implementation() != result.implementation()) {
			vec_repurpose(
				vec.communicator(),
				vec.type(),
				vec.local_size(),
				vec.size(),
				&result.implementation()
			);

			check_error( VecCopy(vec.implementation(), result.implementation()) );
		}	

		VecAbs(result.implementation());
		result.set_initialized(true);
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
	
	void PetscBackend::apply_binary(PetscVector &result, const PetscVector &left, const EMultiplies &, const PetscVector &right) {
		allocate_apply_vec(result, left, right);
		check_error(VecPointwiseMult(result.implementation(), left.implementation(), right.implementation()));
		result.set_initialized(true);
	}
	
	void PetscBackend::apply_binary(PetscVector &result, const PetscVector &left, const Divides &, const PetscVector &right) {
		allocate_apply_vec(result, left, right);
		check_error(VecPointwiseDivide(result.implementation(), left.implementation(), right.implementation()));
		result.set_initialized(true);
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
		if(vec.implementation() != result.implementation()) {
			vec_repurpose(vec.communicator(), vec.type(), vec.local_size(), vec.size(), &result.implementation());
			check_error( VecCopy(vec.implementation(), result.implementation()) );
		}	

		if(reciprocal.numerator() == 1) {
			check_error(VecReciprocal(result.implementation()) );
		} else {
			Vec multiplier;
			check_error( VecDuplicate(result.implementation(), &multiplier) );
			check_error( VecSet(multiplier, reciprocal.numerator()) );
			check_error( VecPointwiseDivide(result.implementation(), multiplier, result.implementation()) );
			check_error( VecDestroy(&multiplier) );
		}

		result.set_initialized(true);
	}
	
	Scalar PetscBackend::norm2(const PetscVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_2, &val);
		return val;
	}
	
	Scalar PetscBackend::norm2(const Matrix &m) {
		PetscReal val;
		MatNorm(m.implementation(), NORM_FROBENIUS, &val);
		return val;
	}
	
	Scalar PetscBackend::norm1(const PetscVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_1, &val);
		return val;
	}
	
	Scalar PetscBackend::norm_infty(const PetscVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_INFINITY, &val);
		return val;
	}
	
	Scalar PetscBackend::reduce(const PetscVector &vec, const Plus &) {
		Scalar result = 0;
		VecSum(vec.implementation(), &result);
		return result;
	}
	
	Scalar PetscBackend::reduce(const PetscMatrix &mat, const Plus &op) {
		PetscVector row_sum;
		//FIXME create compatible vector deduced from matrix type
		vec_create_parallel(mat.communicator(), mat.local_size().get(0), mat.size().get(0), &row_sum.implementation());
		MatGetRowSum(mat.implementation(), row_sum.implementation());
		return reduce(row_sum, op);
	}
	
	Scalar PetscBackend::reduce(const PetscVector &v, const Min &)
	{
		Scalar x = 0.0;
		VecMin(v.implementation(), nullptr, &x);
		return x;
	}

	template<class Operation>
	inline static Scalar generic_local_reduce(const PetscMatrix &m, const Scalar &init_value, const Operation &op)
	{
		Scalar x = init_value;
		const Scalar * values;
		const PetscInt * cols;
		
		PetscInt r_begin, r_end;
		PetscInt n_values = 0;

		PetscInt local_r, local_c;
		MatGetLocalSize(m.implementation(), &local_r, &local_c);
		MatGetOwnershipRange(m.implementation(), &r_begin, &r_end);
		
		for(PetscInt row = r_begin; row < r_end; ++row) {	

			MatGetRow(m.implementation(), row, &n_values, &cols, &values);

			if(n_values < local_c) {
				x = op.template apply<Scalar>(x, 0.);
			}

			for(PetscInt i = 0; i < n_values; ++i) {
				x = op.template apply<Scalar>(x, values[i]);
			}

			MatRestoreRow(m.implementation(), row, &n_values, &cols, &values);
		}

		return x;
	}

	Scalar PetscBackend::reduce(const PetscMatrix &m, const Min &op)
	{
		Scalar x = std::numeric_limits<Scalar>::max();
		MPI_Comm comm = m.communicator();
		
		int size = 0;
		MPI_Comm_size(comm, &size);
		if(size == 1) {
			PetscInt grows, gcols, lrows, lcols;
			MatGetSize(m.implementation(), &grows, &gcols);
			MatGetLocalSize(m.implementation(), &lrows, &lcols);
			
			Vec v;
			//FIXME create compatible vector deduced from matrix type
			vec_create_parallel(comm, lrows, grows, &v);
			MatGetRowMin(m.implementation(), v, nullptr);
			VecMin(v, nullptr, &x);

			VecDestroy(&v);
			return x;
		} else {
			x = generic_local_reduce(m, x, op);
			MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_MIN, comm);
			return x;
		}
	}

	Scalar PetscBackend::reduce(const PetscVector &v, const Max &)
	{
		Scalar x = -std::numeric_limits<Scalar>::max();
		VecMax(v.implementation(), nullptr, &x);
		return x;
	}

	Scalar PetscBackend::reduce(const PetscMatrix &m, const Max &op)
	{
		Scalar x = -std::numeric_limits<Scalar>::max();
		MPI_Comm comm = m.communicator();
		
		int size = 0;
		MPI_Comm_size(comm, &size);
		if(size == 1) {
			PetscInt grows, gcols, lrows, lcols;
			MatGetSize(m.implementation(), &grows, &gcols);
			MatGetLocalSize(m.implementation(), &lrows, &lcols);
			
			Vec v;
			//FIXME create compatible vector deduced from matrix type
			vec_create_parallel(comm, lrows, grows, &v);
			MatGetRowMax(m.implementation(), v, nullptr);
			VecMax(v, nullptr, &x);

			VecDestroy(&v);
			return x;
		} else {
			x = generic_local_reduce(m, x, op);
			MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_MAX, comm);
			return x;
		}
	}

	// get diagonal of matrix as vector
	void PetscBackend::diag(PetscVector &vec, const PetscMatrix &mat)
	{
		using std::max;
			
		auto gs = mat.size();
		auto ls = mat.local_size();
		
		PetscInt len_global = max(gs.get(0), gs.get(1));
		PetscInt len_local  = max(ls.get(0), ls.get(1));

		//FIXME create compatible vector deduced from matrix type
		vec_create_parallel(mat.communicator(), len_local, len_global, &vec.implementation());
		check_error( MatGetDiagonal(mat.implementation(), vec.implementation()) );

		vec.set_initialized(true);
	}
	
	void PetscBackend::diag(PetscSparseMatrix &mat, const PetscVector &vec)
	{
		MPI_Comm comm = vec.communicator();		
		PetscInt local_size  = vec.local_size();
		PetscInt global_size = vec.size();

		check_error( MatDestroy(&mat.implementation()) );
		
		sparse_mat_create_parallel(
			comm,
			local_size,
			local_size,
			global_size,
			global_size,
			1,
			1,
			&mat.implementation());
		
		check_error( MatDiagonalSet( mat.implementation(), vec.implementation(), INSERT_VALUES) );
	}

	void PetscBackend::diag(PetscMatrix &out, const PetscMatrix &in)
	{
		PetscVector vec;
		diag(vec, in);
		diag(out, vec);
	}

	void PetscBackend::diag(PetscMatrix &mat, const PetscVector &vec)
	{
		MPI_Comm comm = vec.communicator();
		PetscInt local_size  = vec.local_size();
		PetscInt global_size = vec.size();
		
		check_error( MatDestroy(&mat.implementation()) );

		dense_mat_create_parallel(
			comm,
			local_size,
			local_size,
			global_size,
			global_size,
			&mat.implementation()
			);
		
		check_error( MatDiagonalSet( mat.implementation(), vec.implementation(), INSERT_VALUES) );
	}

	void PetscBackend::mat_diag_shift(PetscMatrix &left, const Scalar diag_factor)
	{
		check_error(MatShift(left.implementation(), diag_factor));
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
		VecScale(result.implementation(), factor);
	}
	
	void PetscBackend::apply_binary(Matrix &result, const Scalar factor, const Multiplies &, const Matrix &mat) {
		result = mat;
		MatScale(result.implementation(), factor);
	}

	void PetscBackend::scale(Vector &result, const Scalar scale_factor)
	{
		VecScale(result.implementation(), scale_factor);
	}

	void PetscBackend::scale(Matrix &result, const Scalar scale_factor)
	{
		MatScale(result.implementation(), scale_factor);
	}

	Scalar PetscBackend::dot(const Vector &left, const Vector &right) {
		Scalar result;
		VecDot(left.implementation(), right.implementation(), &result);
		return result;
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
			MatDestroy(&result.implementation());
		}
		
		check_error( MatMatMatMult(A.implementation(), B.implementation(), C.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
	}

	bool PetscBackend::is_nan_or_inf(const PetscVector &X)
	{
		PetscInt m; 
		const Scalar *x;
		VecGetLocalSize(X.implementation(), &m);
		VecGetArrayRead(X.implementation(), &x);

		int has_nan = 0; 

		for (PetscInt i = 0; i < m; i++) {
			has_nan = PetscIsInfOrNanScalar(x[i]); 
			if(has_nan == 1)
				break;
		}

		VecRestoreArrayRead(X.implementation(), &x);
		MPI_Comm comm = PetscObjectComm((PetscObject) X.implementation());
		MPI_Allreduce(MPI_IN_PLACE, &has_nan, 1, MPI_INT, MPI_MAX, comm);
		return has_nan > 0; 
	}

	bool PetscBackend::is_nan_or_inf(const PetscMatrix &m)
	{
		int has_nan = 0;
		const Scalar * values;
		const PetscInt * cols;
		
		PetscInt r_begin, r_end;
		PetscInt n_values = 0;

		PetscInt local_r, local_c;
		MatGetLocalSize(m.implementation(), &local_r, &local_c);
		MatGetOwnershipRange(m.implementation(), &r_begin, &r_end);
		
		for(PetscInt row = r_begin; row < r_end; ++row) {	

			MatGetRow(m.implementation(), row, &n_values, &cols, &values);

			for(PetscInt i = 0; i < n_values; ++i) {
				has_nan = (PetscIsInfOrNanScalar(values[i])? 1 : 0);
				if(has_nan) break;
			}

			MatRestoreRow(m.implementation(), row, &n_values, &cols, &values);
			if(has_nan) break;
		}

		MPI_Allreduce(MPI_IN_PLACE, &has_nan, 1, MPI_INT, MPI_MAX, m.communicator());
		return has_nan > 0;
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
		check_error( VecAXPY(y.implementation(), alpha, x.implementation()) );
		y.set_initialized(true);
	}

	void PetscBackend::axpy(PetscMatrix &y, const Scalar alpha, const PetscMatrix &x)
	{
		check_error( MatAXPY(y.implementation(), alpha, x.implementation(), DIFFERENT_NONZERO_PATTERN) );
	}

	void PetscBackend::axpby(Vector &y, const Scalar alpha, const Vector &x, const Scalar &beta)
	{
		check_error( VecAXPBY(y.implementation(), alpha, beta, x.implementation()) );
		y.set_initialized(true);
	}
	
	void PetscBackend::diag_scale_right(Matrix &result, const Matrix &m, const Vector &diag)
	{
		assign(result, m);
		check_error( MatDiagonalScale(result.implementation(), nullptr, diag.implementation()) );
	}

	void PetscBackend::diag_scale_left(Matrix &result, const Vector &diag, const Matrix &m)
	{
		assign(result, m);
		check_error( MatDiagonalScale(result.implementation(), diag.implementation(), nullptr) );
	}

	void PetscBackend::diag_scale_left(Vector &result, const Vector &diag, const Vector &m)
	{
		if(result.implementation() != diag.implementation() && result.implementation() != m.implementation()) {
			vec_repurpose(diag.communicator(), diag.type(), diag.local_size(), diag.size(), &result.implementation());
		}

		check_error( VecPointwiseMult(diag.implementation(), m.implementation(), result.implementation()) );
		result.set_initialized(true);
	}

	Scalar PetscBackend::trace(const Matrix &mat)
	{
		Scalar ret;
		check_error( MatGetTrace(mat.implementation(), &ret) );
		return ret;
	}

	template<class Operation>
	inline static void generic_col_reduce(
		PetscVector &result,
		const PetscMatrix &mat, 
		const Scalar &init_value, 
		const Operation &op
		)
	{
		
		const Scalar * values;
		const PetscInt * cols;
		
		PetscInt r_begin, r_end;
		PetscInt n_values = 0;

		PetscInt global_r, global_c, local_r, local_c;
		MatGetSize(mat.implementation(), &global_r, &global_c);
		MatGetLocalSize(mat.implementation(), &local_r, &local_c);

		VecDestroy(&result.implementation());
		PetscBackend::vec_create_parallel(mat.communicator(), local_r, global_r, &result.implementation());
		MatGetOwnershipRange(mat.implementation(), &r_begin, &r_end);
		
		PetscBackend::write_lock(result);

		for(PetscInt row = r_begin; row < r_end; ++row) {	
			MatGetRow(mat.implementation(), row, &n_values, &cols, &values);

			Scalar x = init_value;
			for(PetscInt i = 0; i < n_values; ++i) {
				x = op.template apply<Scalar>(x, values[i]);
			}

			if(n_values < local_c) {
				x = op.template apply<Scalar>(x, 0.);
			}

			MatRestoreRow(mat.implementation(), row, &n_values, &cols, &values);
			VecSetValues(result.implementation(), 1, &row, &x, INSERT_VALUES);
		}

		PetscBackend::write_unlock(result);
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Plus &, const int dim)
	{
		vec_repurpose(
			mat.communicator(),
			parallel_vector_type(),
			mat.local_size().get(0),
			mat.size().get(0),
			&result.implementation());

		if(dim == 1) {
			MatGetRowSum(mat.implementation(), result.implementation());
			result.set_initialized(true);
		} else {
			assert(false && "not available in pestsc");
			result.set_initialized(false);
		}
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Min &op, const int dim)
	{
		PetscInt grows, gcols;
		MatGetSize(mat.implementation(), &grows, &gcols);

		int comm_size = 0;
		MPI_Comm_size(mat.communicator(), &comm_size);
		bool is_serial = comm_size == 1;

		if (dim == 1) {
			if(is_serial) {
				if(!result.is_null()) {
					VecDestroy(&result.implementation());
				}

				vec_create_parallel(mat.communicator(), PETSC_DECIDE, grows, &result.implementation());
				MatGetRowMin(mat.implementation(), result.implementation(), nullptr);
			} else {
				generic_col_reduce(
					result,
					mat, 
					std::numeric_limits<Scalar>::max(), 
					op
					);

				result.set_initialized(true);
			}
		} else {
			assert(false && "not available in petsc");
			result.set_initialized(false);
		}
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Max &op, const int dim)
	{
		PetscInt grows, gcols;
		MatGetSize(mat.implementation(), &grows, &gcols);

		int comm_size = 0;
		MPI_Comm_size(mat.communicator(), &comm_size);
		bool is_serial = comm_size == 1;

		if (dim == 1) {
			if(is_serial) {
				check_error( VecDestroy(&result.implementation()) );
				vec_create_parallel(mat.communicator(), PETSC_DECIDE, grows, &result.implementation());
				check_error( MatGetRowMax(mat.implementation(), result.implementation(), nullptr) );
			} else {
				generic_col_reduce(
					result,
					mat, 
					-std::numeric_limits<Scalar>::max(), 
					op);

				result.set_initialized(true);
			}
		} else {
			assert(false && "not available in petsc");
			result.set_initialized(false);
		}

	}

	void PetscBackend::inverse(Matrix &result, const Matrix &mat)
	{
		PetscMatrix I;

		Size s;
		size(mat, s);

		build(I, s, Identity());
		build(result, s, Zeros());

		IS isr, isc;
		MatFactorInfo info;

		Matrix L;
		assign(L, mat);

		check_error( MatGetOrdering(L.implementation(), MATORDERINGNATURAL, &isr, &isc) );
		check_error( MatLUFactor( L.implementation(), isr, isc, &info ) ); 
		check_error( MatMatSolve(L.implementation(), I.implementation(), result.implementation()) );
		check_error( ISDestroy(&isr) );
		check_error( ISDestroy(&isc) );
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
			PetscObjectSetName((PetscObject)m.implementation(), args.name.c_str());
		}
	}

	void PetscBackend::apply_args(const PetscArgs &args, Vector &v)
	{
		if(!args.name.empty()) {
			PetscObjectSetName((PetscObject)v.implementation(), args.name.c_str());
		}
	}

	void PetscBackend::build_ghosts(
		const PetscInt &local_size,
		const PetscInt &global_size,
		const std::vector<PetscInt> &index,
		PetscVector &vec)
	{
		VecDestroy(&vec.implementation());	
		vec.set_initialized(false);

		check_error(
			VecCreateGhost(
				default_communicator(), 
				local_size,
				global_size,
				static_cast<PetscInt>(index.size()),
				&index[0],
				&vec.implementation())
		);

		vec.init_ghost_index(index);

		VecZeroEntries(vec.implementation());
	}

	void PetscBackend::update_ghosts(PetscVector &vec)
	{
		vec.update_ghosts();
	}

}

