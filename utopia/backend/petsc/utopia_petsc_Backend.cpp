
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
					PETScError::Check( MatGetType(right, &common_type) );
					PETScError::Check( MatConvert(left, common_type, MAT_INITIAL_MATRIX, &left_) );
					right_ = right;
					must_destroy_left_  = true;
			
				} else {
					MatCreate(comm, &right_);
					PETScError::Check( MatGetType(left, &common_type) );
					PETScError::Check( MatConvert(right, common_type, MAT_INITIAL_MATRIX, &right_) );
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

	void PetscBackend::resize(PETScMatrix &mat, const Size &s)
	{
		Mat &m = mat.implementation();
		PetscInt r, c;
		MatGetSize(mat.implementation(), &r, &c);
		if(r == s.get(0) && c == s.get(1)) {
			return;
		}

		MPI_Comm comm = PetscObjectComm((PetscObject) mat.implementation());
		MatType type;
		PETScError::Check( MatGetType(m, &type) );
		MatDestroy(&m);

		MatCreate(comm, &m);
		MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1));
		MatSetType(m, type);
	}

	void PetscBackend::resize(PETScMatrix &mat, const Size &local_s, const Size &global_s)
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
		PETScError::Check( MatGetType(m, &type) );
		MatDestroy(&m);

		MatCreate(comm, &m);
		MatSetSizes(m, local_s.get(0), local_s.get(1), global_s.get(0), global_s.get(1));
		MatSetType(m, type);
	}

	void PetscBackend::resize(PETScVector &vec, const Size &s)
	{
		Vec &v = vec.implementation();
		PetscInt n;

		if(vec.isInitialized()) {
			VecGetSize(v,  &n);
			if(n == s.get(0)) {
				return;
			}
		}

		MPI_Comm comm = PetscObjectComm((PetscObject) vec.implementation());
		VecDestroy(&v);	
		VecCreateMPI(comm, PETSC_DECIDE, s.get(0), &v);
	}


	void PetscBackend::resize(PETScVector &vec, const Size &s_local, const Size &s_global)
	{
		Vec &v = vec.implementation();
		PetscInt n, N;

		if(vec.isInitialized()) {
			VecGetSize(v, &N);
			VecGetLocalSize(v,  &n);
			if(n == s_local.get(0) && N == s_global.get(0)) {
				return;
			}
		}

		MPI_Comm comm = PetscObjectComm((PetscObject) vec.implementation());
		VecDestroy(&v);	
		VecCreateMPI(comm, s_local.get(0), s_global.get(0), &v);
	}

	void PetscBackend::assign_transposed(PETScMatrix &left, const PETScMatrix &right) {
		if(&left != &right) {
			MatDestroy(&left.implementation());
		} else {
			assert(false);
		}

		MatTranspose(right.implementation(), MAT_INITIAL_MATRIX, &left.implementation());
	}
	
	void PetscBackend::clear(PETScMatrix &mat)
	{
		MatDestroy(&mat.implementation());
		MatCreate(mat.communicator(), &mat.implementation());
	}
	
	void PetscBackend::convert(Vec vec, PETScVector &wrapper)
	{
		MPI_Comm comm = PetscObjectComm((PetscObject) vec);

		PetscInt n, local_n;
		VecGetSize(vec, &n);
		VecGetLocalSize(vec, &local_n);
		VecDestroy(&wrapper.implementation());
		VecCreateMPI(comm, local_n, n, &wrapper.implementation());
		VecCopy(vec, wrapper.implementation());
	}
	
	void PetscBackend::convert(Mat mat, PETScMatrix &wrapper)
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
	
	void PetscBackend::convert(Mat mat, PETScSparseMatrix &wrapper)
	{
		MatDestroy(&wrapper.implementation());
		MatDuplicate(mat, MAT_COPY_VALUES, &wrapper.implementation());
	}
	
	void PetscBackend::convert(const PETScVector &wrapper, Vec vec)
	{
		VecCopy(wrapper.implementation(), vec);
	}
	
	void PetscBackend::convert(const PETScMatrix &wrapper, Mat mat)
	{
		MatCopy(wrapper.implementation(), mat, DIFFERENT_NONZERO_PATTERN);
	}
	
	void PetscBackend::convert(const PETScSparseMatrix &wrapper, Mat mat)
	{
		MatCopy(wrapper.implementation(), mat, DIFFERENT_NONZERO_PATTERN);
	}

	void PetscBackend::wrap(Mat mat, PETScSparseMatrix &wrapper)
	{
		wrapper.wrap(mat);
	}

	void PetscBackend::wrap(Vec vec, PETScVector &wrapper)
	{
		assert(false && "TODO");
		// wrapper.wrap(vec);
	}

	Range PetscBackend::range(const PETScVector &v) 
	{
		PetscInt rbegin, rend;
		VecGetOwnershipRange(v.implementation(), &rbegin, &rend);
		return Range(rbegin, rend);
	}
	
	Range PetscBackend::row_range(const PETScMatrix &m) 
	{
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		assert(Range(rbegin, rend).valid());
		return Range(rbegin, rend);
	}
	
	Range PetscBackend::col_range(const PETScMatrix &m) 
	{
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		return Range(0, gcols);
	}
	
	void PetscBackend::size(const PETScMatrix &m, Size &size) 
	{
		PetscInt grows, gcols;
		size.set_dims(2);
		MatGetSize(m.implementation(), &grows, &gcols);
		size.set(0, grows);
		size.set(1, gcols);
	}
	
	void PetscBackend::size(const PETScVector &m, Size &size) 
	{
		PetscInt n;
		size.set_dims(1);
		VecGetSize(m.implementation(), &n);
		size.set(0, n);
	}
	
	void PetscBackend::local_size(const PETScVector &m, Size &size) 
	{
		PetscInt n;
		size.set_dims(1);
		VecGetLocalSize(m.implementation(), &n);
		size.set(0, n);
	}

	void PetscBackend::local_size(const PETScMatrix &mat, Size &size)
	{
		PetscInt n, m;
		size.set_dims(2);
		MatGetLocalSize(mat.implementation(), &n, &m);
		size.set(0, n);
		size.set(1, m);
	}

	void PetscBackend::select(PETScVector &left,
				const PETScVector &right,
	      		const std::vector<PetscInt> &index)
	{

		Vec l = left.implementation();
		Vec r = right.implementation();

		MPI_Comm comm = PetscObjectComm((PetscObject)r);

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


	void PetscBackend::select(PETScMatrix &left,
				const PETScMatrix &right,
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

	void PetscBackend::select_aux(PETScMatrix &left,
				const PETScMatrix &right,
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
		const PETScMatrix &right,
		PETScMatrix &result)
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
		const PETScMatrix &right,
		PETScMatrix &result)
	{

		Mat &l = result.implementation();
		const Mat r = right.implementation();


		std::stringstream ss;

		for(auto r : remote_rows) {
			ss << r;
		}

		IS isrow;
		PetscErrorCode ierr = ISCreateGeneral(PETSC_COMM_WORLD, remote_rows.size(), &remote_rows[0], PETSC_USE_POINTER, &isrow);

		IS iscol;
		ierr = ISCreateGeneral(PETSC_COMM_WORLD, remote_cols.size(), &remote_cols[0], PETSC_USE_POINTER, &iscol);

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
		const PETScMatrix &right,
		PETScMatrix &result)
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
		PETScMatrix &left,
		const PETScMatrix &right,
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
	bool PetscBackend::read(const std::string &path, PETScMatrix &Mat_A)
	{
			
		Mat &A = Mat_A.implementation();
		MatDestroy(&A);
		PetscViewer fd;
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_READ, &fd);
		MatCreate(PETSC_COMM_WORLD,&A);
		bool status;
		status =  PETScError::Check( MatLoad(A,fd) );
		PetscViewerDestroy(&fd);
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
	bool PetscBackend::write(const std::string &path, const PETScMatrix &Mat_A)
	{
		const Mat &A = Mat_A.implementation();

		const bool is_matlab = is_matlab_file(path);
				
		if(is_matlab) {
			PetscObjectSetName((PetscObject)A, "matrix");

			PetscErrorCode ierr;
			PetscViewer fd;
			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,path.c_str(), &fd); //CHKERRV(ierr);
			ierr = PetscViewerPushFormat(fd,PETSC_VIEWER_ASCII_MATLAB); //CHKERRV(ierr);
			ierr = MatView(A, fd); //CHKERRV(ierr);
			PetscViewerDestroy(&fd);
			return PETScError::Check(ierr);
		} else {
			PetscViewer fd;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_WRITE, &fd);
			bool status;
			status =  PETScError::Check( MatView(A,fd));
			PetscViewerDestroy(&fd);
			return status;
		}
	}
	
	// write vector
	bool PetscBackend::write(const std::string &path, const PETScVector &Vec_A)
	{
		const Vec &A = Vec_A.implementation();

		bool is_matlab = is_matlab_file(path);
		if(is_matlab) {
			PetscObjectSetName((PetscObject)A, "vector");

			PetscErrorCode ierr;
			PetscViewer fd;
			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,path.c_str(), &fd); //CHKERRV(ierr);
			ierr = PetscViewerPushFormat(fd,PETSC_VIEWER_ASCII_MATLAB); //CHKERRV(ierr);
			ierr = VecView(A, fd); //CHKERRV(ierr);
			PetscViewerDestroy(&fd);
			return PETScError::Check(ierr);
		} else {
			PetscViewer fd;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_WRITE, &fd);
			bool status;
			status =  PETScError::Check( VecView(A,fd));
			PetscViewerDestroy(&fd);
			return status;
		}
	}
	


	// monitor for cyrill  - maybe u need to adjust something here ...
	void PetscBackend::monitor(const long & it, PETScMatrix &Mat_A)
	{
		PetscViewer viewer_hessian = nullptr;
		if(it==0)
		{
          	// log stifness
          	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "log_hessian.m" ,&viewer_hessian);  
          	PetscViewerPushFormat(viewer_hessian,PETSC_VIEWER_ASCII_MATLAB); 
        }   

        //FIXME why are you duplicating and then copying?
    	Mat A; 
    	MatDuplicate(Mat_A.implementation(), MAT_COPY_VALUES,  &A); 
    	MatCopy(Mat_A.implementation(), A, SAME_NONZERO_PATTERN); 

	    PetscObjectSetName((PetscObject)A, "H");
	    MatView(A, viewer_hessian); 
		//PetscViewerDestroy(&viewer_hessian);
	}
	
	// monitor for cyrill  - maybe u need to adjust something here ...
	void PetscBackend::monitor(const long & it, PETScVector &Vec_A)
	{
		PetscViewer viewer_iterates = nullptr;
		if(it==0)
		{
          	// log iterates
          	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "log_iterate.m", &viewer_iterates);  
          	PetscViewerPushFormat(viewer_iterates,PETSC_VIEWER_ASCII_MATLAB); 
        }

        //FIXME why are you duplicating and then copying?
        Vec iterates; 	        
        VecDuplicate(Vec_A.implementation(), &iterates);
  		VecCopy(Vec_A.implementation(),iterates);
  		PetscObjectSetName((PetscObject)iterates,"it");
  		VecView(iterates, viewer_iterates); 

      	//PetscViewerDestroy(&viewer_iterates);
	}

	Scalar PetscBackend::get_global_nnz(PETScMatrix &Mat_A)
	{	
		MatInfo        info;
		MatGetInfo(Mat_A.implementation(), MAT_GLOBAL_SUM, &info); 
		return info.nz_used; 
	}

	Scalar PetscBackend::get_local_nnz(PETScMatrix &Mat_A)
	{	
		MatInfo        info;
		MatGetInfo(Mat_A.implementation(), MAT_LOCAL, &info); 
		return info.nz_used; 
	}

	void PetscBackend::set_zero_rows(PETScMatrix &Mat_A, const std::vector<int> &index)
	{
		PETScError::Check(MatZeroRows(Mat_A.implementation(), index.size(), &index[0], 1.0, NULL, NULL));
	}

	void PetscBackend::apply_BC_to_system(PETScMatrix & A, PETScVector& x, PETScVector& rhs, const std::vector<int> &index)
	{
		PETScError::Check(MatZeroRows(A.implementation(), index.size(), &index[0], 1.0, x.implementation(), rhs.implementation())); 
	}

	// read vector
	bool PetscBackend::read(const std::string &path, PETScVector &Vec_A)
	{
		MPI_Comm comm = Vec_A.communicator();

		Vec &A = Vec_A.implementation();
		VecDestroy(&A);

		PetscViewer fd;
		PetscViewerBinaryOpen(comm, path.c_str(), FILE_MODE_READ, &fd);
		VecCreate(comm, &A);
		bool status;
		status =  PETScError::Check( VecLoad(A, fd));
		PetscViewerDestroy(&fd);
		return status;
	}
	
	void PetscBackend::assign_from_range(
		PETScVector &left, 
		const PETScVector &right, 
		const Range &global_row_range,
		const Range &) 
	{
		assert(!global_row_range.empty());
		
		Vec &li = left.implementation();
		//const Vec &ri = left.implementation();
		
		// get range of vector
		Range rr = range(right).intersect(global_row_range);
		
		// create vector
		//VecCreate(left.communicator(), &li);
		// TODO Check VecCreate on already exisiting vec
		VecSetType(li, VECMPI);
		
		// set vector size
		VecSetSizes(li, PETSC_DECIDE, rr.extent());
		
		
		
		PetscInt r = 0;
		for (PetscInt rIt = rr.begin(); rIt < rr.end(); ++rIt) {
			r = rIt - global_row_range.begin();
			set(left, r, get(right, rIt));
		}
		
		VecAssemblyBegin(li);
		VecAssemblyEnd(li);
	}
	
	void PetscBackend::gemm(
		PETScMatrix &result,
		const Scalar beta,
		const Scalar alpha,
		bool transpose_left,
		const PETScMatrix &left,
		bool transpose_right,
		const PETScMatrix &right) 
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
			ok = PETScError::Check(MatTransposeMatMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()));
		} else if(!transpose_left && transpose_right) {
			ok = PETScError::Check(MatMatTransposeMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()));
		} else if(!transpose_left && !transpose_right) {
			ok = PETScError::Check(MatMatMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()));
		} else {
			assert(transpose_left && transpose_right);
			PETScMatrix temp;

			if(!PETScError::Check( MatMatMult(r, l, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation()) )) {
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
		const PETScMatrix &A, 
		const Vector &x)
	{
		assert(fabs(beta) < 1e-16);
		assert(fabs(alpha - 1.) < 1e-16);

		const Mat &A_im = A.implementation();
		Vec &y_im = y.implementation();
		const Vec &x_im = x.implementation();
		
		 Size gs, ls;
		 size(A, gs);
		 local_size(A, ls);

		 VecType type;
		 VecGetType(x_im, &type);
		 VecSetType(y_im, type);
		 
		 if(transpose_A) {
		 	VecSetSizes(y_im, ls.get(1), gs.get(1));
		 	PETScError::Check( MatMultTranspose(A_im, x_im, y_im) );
		 } else {
		 	VecSetSizes(y_im, ls.get(0), gs.get(0));
		 	PETScError::Check( MatMult(A_im, x_im, y_im) );
		 }
	}
		
	void PetscBackend::build(PETScMatrix &m, const Size &size, const Identity &) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		//FIXME use this: MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],Scalar diag,Vec x,Vec b)
		MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), NULL,
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
	}
	
	void PetscBackend::build(PETScSparseMatrix &m, const Size &size, const Identity &) {
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, rows, cols, 1, PETSC_NULL,
					 1 /*Only because otherwise petsc crashes*/, PETSC_NULL, &m.implementation());
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE);
		MatSetOption(m.implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
		MatSetOption(m.implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);

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
	}
	
	void PetscBackend::build(PETScMatrix &m, const Size &size, const LocalIdentity &) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		//FIXME use this: MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],Scalar diag,Vec x,Vec b)
		MatCreateDense(comm, size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, NULL,
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
	}
	
	void PetscBackend::build(PETScSparseMatrix &m, const Size &size, const LocalIdentity &) {
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, rows, cols, PETSC_DETERMINE, PETSC_DETERMINE, 1, PETSC_NULL,
					 1 /*Only because otherwise petsc crashes*/, PETSC_NULL, &m.implementation());
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE);
		MatSetOption(m.implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
		MatSetOption(m.implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);

		unsigned long send_buffer = cols;
		unsigned long receive_buffer = 0;
		
		MPI_Exscan(&send_buffer, &receive_buffer, 1, MPI_UNSIGNED_LONG ,
				   MPI_SUM, m.communicator());
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt extent = PetscMin(rend-rbegin, PetscMin(rows, cols));

		PETScError::Check( MatZeroEntries(m.implementation()) );
		
		for (PetscInt i = 0; i < extent; ++i) {
			const PetscInt row_offset = rbegin+i;
			const PetscInt col_offset = receive_buffer+i;
			const Scalar val = 1;
			MatSetValues(m.implementation(), 1, &row_offset, 1, &col_offset, &val, INSERT_VALUES);
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PetscBackend::build(PETScSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz) {
		
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, rows, cols,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to local entries*/,  PETSC_NULL,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to remote entries*/, PETSC_NULL,
					 &m.implementation());

		// MatSetOption(m.implementation(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE);
		MatSetOption(m.implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
		MatSetOption(m.implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);

		PETScError::Check( MatZeroEntries(m.implementation()) );
	}
	
	void PetscBackend::build(PETScSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz) {
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		Mat &mat = m.implementation();

		MatDestroy(&mat);
		MatCreate(comm, &mat);
		MatSetSizes(mat, rows, cols, PETSC_DETERMINE, PETSC_DETERMINE);
		
		PETScError::Check( MatSetType(mat, MATAIJ) );
		PETScError::Check( MatSeqAIJSetPreallocation(mat, PetscMax(nnz.nnz(), 1), PETSC_NULL) );
		PETScError::Check( MatMPIAIJSetPreallocation(mat, PetscMax(nnz.nnz(), 1), PETSC_NULL, PetscMax(nnz.nnz(), 1), PETSC_NULL) ); 
					
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE);
		MatSetOption(m.implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
		MatSetOption(m.implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);

		PETScError::Check( MatZeroEntries(m.implementation()) );
	}
	
	void PetscBackend::build(PETScSparseMatrix &m, const Size &size, const LocalRowNNZ<PetscInt> &nnz)
	{
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, rows, PETSC_DECIDE, PETSC_DETERMINE, cols,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to local entries*/, PETSC_NULL,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to remote entries*/, PETSC_NULL,
					 &m.implementation());
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE);
		MatSetOption(m.implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
		MatSetOption(m.implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);

		PETScError::Check( MatZeroEntries(m.implementation()) );
	}
	
	/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
	void PetscBackend::build(PETScMatrix  &m, const Size &size, const LocalNNZ<PetscInt> & /*nnz */)
	{
		build(m, size, LocalValues<Scalar>(0));
	}
	
	/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
	void PetscBackend::build(PETScMatrix  &m, const Size &size, const NNZ<PetscInt> &/*nnz*/)
	{
		build(m, size, Zeros());
	}
	
	void PetscBackend::build(PETScMatrix &m, const Size &size, const Zeros &)
	{
		build(m, size, Values<Scalar>(0));
	}
	
	void PetscBackend::build(PETScVector &v, const Size &size, const Zeros &)
	{
		build(v, size, Values<Scalar>(0));
	}
	
	void PetscBackend::build(PETScMatrix &m, const Size &size, const LocalZeros &)
	{
		build(m, size, LocalValues<Scalar>(0));
	}
	
	void PetscBackend::build(PETScVector &v, const Size &size, const LocalZeros &)
	{
		build(v, size, LocalValues<Scalar>(0));
	}
	
	void PetscBackend::build(PETScMatrix &m, const Size &size, const Values<Scalar> &values) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), NULL,
					   &m.implementation());
		
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		PETScError::Check( MatZeroEntries(m.implementation()) );
		
		
		
		const Scalar v = values.value();
		for (PetscInt i = rbegin; i < rend; ++i) {
			for (PetscInt j = 0; j < gcols; ++j) {
				MatSetValue(m.implementation(), i, j, v, INSERT_VALUES);
			}
		}
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PetscBackend::build(PETScVector &v, const Size &in_local_size, const Size &&in_global_size, const Values<Scalar> &values)
	{
		Size v_local_size;
		if(v.isInitialized()) {
			local_size(v, v_local_size);
		}	

		if(!v.isInitialized() || v_local_size.get(0) != in_local_size.get(0)) {
			VecDestroy(&v.implementation());
			VecCreateMPI(v.communicator(), in_local_size.get(0), in_global_size.get(0), &v.implementation());
		}
 
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());
	}
	
	void PetscBackend::build(PETScVector &v, const Size &size, const Values<Scalar> &values) {

		VecDestroy(&v.implementation());
		VecCreateMPI(v.communicator(), PETSC_DECIDE, size.get(0), &v.implementation());
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());
	}
	
	void PetscBackend::build(PETScMatrix &m, const Size &size, const LocalValues<Scalar> &values) {
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
	}
	
	void PetscBackend::build(PETScVector &v, const Size &size, const LocalValues<Scalar> &values) {
		MPI_Comm comm = v.communicator();
		VecDestroy(&v.implementation());

		VecCreateMPI(comm, size.get(0), PETSC_DETERMINE, &v.implementation());
		
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());
	}
	
	void PetscBackend::add(PETScVector &v, const PetscInt index, Scalar value)
	{
		VecSetValues(v.implementation(), 1, &index, &value, ADD_VALUES);
	}
	
	void PetscBackend::add(PETScMatrix &m, const PetscInt row, const PetscInt col, Scalar value)
	{
		MatSetValues(m.implementation(), 1, &row, 1, &col, &value, ADD_VALUES);
	}
	
	void PetscBackend::set(PETScVector &v, const PetscInt index, Scalar value) {
		VecSetValues(v.implementation(), 1, &index, &value, INSERT_VALUES);
	}
	
	void PetscBackend::set(PETScVector &v, const std::vector<PetscInt> &indices, const std::vector<Scalar> &values) {
		assert(indices.size() == values.size());
		VecSetValues(v.implementation(), indices.size(), &indices[0], &values[0], INSERT_VALUES);
	}
	
	
	void PetscBackend::set(PETScMatrix &v, const PetscInt row, const PetscInt col, Scalar value) {
		MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
	}
	
	void PetscBackend::set(PETScSparseMatrix &v, const PetscInt row, const PetscInt col, Scalar value) {
		MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
	}
	
	void PetscBackend::write_lock(PETScVector &vec) {}
	
	void PetscBackend::write_unlock(PETScVector &vec) {
		vec.assemblyBegin();
		vec.assemblyEnd();
	}
	
	void PetscBackend::write_lock(const PETScMatrix &) { }
	
	void PetscBackend::write_unlock(const PETScMatrix &mat) {
		MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PetscBackend::write_lock(const PETScSparseMatrix &) { }
	
	void PetscBackend::write_unlock(const PETScSparseMatrix &mat) {
		MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PetscBackend::set(
		PETScMatrix &m,
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
		PETScMatrix &m,
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


	Scalar PetscBackend::get(const PETScVector &v, const PetscInt index) {
		Scalar value;
		VecGetValues(v.implementation(), 1, &index, &value);
		return value;
	}
	
	Scalar PetscBackend::get(const PETScMatrix &v, const PetscInt row, const PetscInt col) {
		Scalar value;
		MatGetValues(v.implementation(), 1, &row, 1, &col, &value);
		return value;
	}
	
	void PetscBackend::apply_binary(PETScVector &result, const PETScMatrix &left, const Multiplies &, const PETScVector &right) {
		PetscInt grows, gcols;
		MatGetSize(left.implementation(), &grows, &gcols);
		
		PetscInt rows, cols;
		MatGetLocalSize(left.implementation(), &rows, &cols);
		
		result.setCommunicator(left.communicator());
		
		VecDestroy(&result.implementation());
		VecCreateMPI(right.communicator(), rows, grows, &result.implementation());
		VecAssemblyBegin(result.implementation());
		VecAssemblyEnd(result.implementation());
		
		MatMult(left.implementation(), right.implementation(), result.implementation());
	}
	
	void PetscBackend::apply_binary(PETScMatrix &result, const PETScMatrix &left, const Multiplies &, const PETScMatrix &right)
	{
		if(&right.implementation() != &result.implementation() || &left.implementation() !=  &result.implementation()) {
			MatDestroy(&result.implementation());
		} else {
			assert(false);
		}
 
		MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
	}

	// void PetscBackend::apply_binary(PETScMatrix &result, const PETScMatrix &left, const Multiplies &, const PETScMatrix &right) {
	// 	Size lsize, rsize;
	// 	size(left, lsize);
	// 	size(right, rsize);
		
	// 	result.initGlobal(lsize.get(0), rsize.get(1));
		
	// 	MatType type;
	// 	MatGetType(right.implementation(), &type);
		
	// 	//            //FIXME: When petsc is fixed for MATMPIDENSE remove branch
	// 	//            if (std::string(MATMPIDENSE) == type) {
	// 	//                Mat sparseleft, sparseright, sparseresult;
	// 	//                MatConvert(left.implementation(), MATMPIAIJ, MAT_INITIAL_MATRIX, &sparseleft);
	// 	//                MatConvert(right.implementation(), MATMPIAIJ, MAT_INITIAL_MATRIX, &sparseright);
	// 	//                if (!PETScError::Check(
	// 	//                        MatMatMult(sparseleft, sparseright, MAT_INITIAL_MATRIX,
	// 	//                                   PETSC_DEFAULT, &sparseresult))) {
	// 	//                    return false;
	// 	//                }
	// 	//
	// 	//                MatConvert(sparseresult, MATMPIDENSE, MAT_INITIAL_MATRIX, &result.implementation());
	// 	//                return true;
	// 	//
	// 	//            }
		
	// 	PETScError::Check( 
	// 		MatMatMult(left.implementation(),
	// 				   right.implementation(),
	// 				   MAT_INITIAL_MATRIX,
	// 				   PETSC_DEFAULT,
	// 				   &result.implementation()));
	// }

	void PetscBackend::mat_mult_add(PETScVector &result, const PETScMatrix &m, const PETScVector &right, const PETScVector &left)
	{
		if (&right.implementation() == &result.implementation() || &left.implementation() ==  &result.implementation()) {
			assert(false);
		}

		PetscInt size, gsize;
		VecGetSize(right.implementation(), &gsize);
		VecGetLocalSize(right.implementation(), &size);

		result.setCommunicator(right.communicator());

		VecDestroy(&result.implementation());
		VecCreateMPI(right.communicator(), size, gsize, &result.implementation());
		VecAssemblyBegin(result.implementation());
		VecAssemblyEnd(result.implementation());

		MatMultAdd(m.implementation(), right.implementation(), left.implementation(), result.implementation());
	}

	void PetscBackend::mat_mult_t_add(PETScVector &result, const PETScMatrix &m, const PETScVector &right, const PETScVector &left)
	{
		if (&right.implementation() == &result.implementation() || &left.implementation() ==  &result.implementation()) {
			assert(false);
		}

		PetscInt size, gsize;
		VecGetSize(right.implementation(), &gsize);
		VecGetLocalSize(right.implementation(), &size);

		result.setCommunicator(right.communicator());

		VecDestroy(&result.implementation());
		VecCreateMPI(right.communicator(), size, gsize, &result.implementation());
		VecAssemblyBegin(result.implementation());
		VecAssemblyEnd(result.implementation());

		MatMultTransposeAdd(m.implementation(), right.implementation(), left.implementation(), result.implementation());
	}

	void PetscBackend::apply_unary(Vector &result, const Abs &, const Vector &vec)
	{
		if(&vec != &result) {
			// bool keep_size = false;
			Size l_s, g_s;

			local_size(vec, l_s);
			size(vec, g_s);

			resize(result, l_s, g_s);
			PETScError::Check( VecCopy(vec.implementation(), result.implementation()) );
			
		}	

		VecAbs(result.implementation());
	}

	void PetscBackend::allocate_apply_vec(PETScVector &result, const PETScVector &left, const PETScVector &right)
	{
		if(&result.implementation() != &left.implementation() && &result.implementation() != &right.implementation()) {
			PetscInt n_wanted; //, n;
			VecGetLocalSize(right.implementation(), &n_wanted);

			// VecGetSize(result.implementation(), &n);

			// if(n != INVALID_INDEX) {
			// 	VecGetLocalSize(result.implementation(), &n);
			// }

			// if(n != n_wanted) {
				PetscInt n_global;
				VecGetSize(right.implementation(), &n_global);
				result.setCommunicator(left.communicator());

				VecDestroy(&result.implementation());
				VecCreateMPI(PetscObjectComm((PetscObject)right.implementation()), n_wanted, n_global, &result.implementation());
			// }
		}
	}
	
	void PetscBackend::apply_binary(PETScVector &result, const PETScVector &left, const EMultiplies &, const PETScVector &right) {
		allocate_apply_vec(result, left, right);
		PETScError::Check(VecPointwiseMult(result.implementation(), left.implementation(), right.implementation()));
	}
	
	void PetscBackend::apply_binary(PETScVector &result, const PETScVector &left, const Divides &, const PETScVector &right) {
		allocate_apply_vec(result, left, right);
		PETScError::Check(VecPointwiseDivide(result.implementation(), left.implementation(), right.implementation()));
	}

	void PetscBackend::apply_binary(PETScVector &result, const PETScVector &left, const Min &op, const PETScVector &right)
	{
		apply_binary_generic(result, left, op, right);
	}

	void PetscBackend::apply_binary(PETScVector &result, const PETScVector &left, const Max &op, const PETScVector &right)
	{
		apply_binary_generic(result, left, op, right);
	}

	// reciprocal
	void PetscBackend::apply_binary(PETScVector &result, const Reciprocal<Scalar> &reciprocal, const PETScVector &vec)
	{
		PetscErrorCode err = 0;

		if(&vec != &result) {
			VecDestroy(&result.implementation());
		 	err = VecDuplicate(vec.implementation(), &result.implementation());
			err = VecCopy(vec.implementation(), result.implementation());
		}	

		if(reciprocal.numerator() == 1)
		{
			err = VecReciprocal(result.implementation());
		} else {
			Vec multiplier;
			
			err =   VecDuplicate(result.implementation(), &multiplier);
			err =   VecSet(multiplier, reciprocal.numerator());
			err =   VecPointwiseDivide(result.implementation(), multiplier, result.implementation());

			VecDestroy(&multiplier);
		}

		PETScError::Check(err);
	}
	
	Scalar PetscBackend::norm2(const PETScVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_2, &val);
		return val;
	}
	
	Scalar PetscBackend::norm2(const Matrix &m) {
		PetscReal val;
		MatNorm(m.implementation(), NORM_FROBENIUS, &val);
		return val;
	}
	
	Scalar PetscBackend::norm1(const PETScVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_1, &val);
		return val;
	}
	
	Scalar PetscBackend::norm_infty(const PETScVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_INFINITY, &val);
		return val;
	}
	
	Scalar PetscBackend::reduce(const PETScVector &vec, const Plus &) {
		//if(!vec.isInitialized()) return 0; // FIXME
		
		Scalar result = 0;
		VecSum(vec.implementation(), &result);
		return result;
	}
	
	Scalar PetscBackend::reduce(const PETScMatrix &mat, const Plus &op) {
		// assert(false);
		PETScVector rowSum; //FIXME initialize
		Size gs, ls;
		size(mat, gs);
		local_size(mat, ls);
		resize(rowSum, ls, gs);

		MatGetRowSum(mat.implementation(), rowSum.implementation());
		return reduce(rowSum, op);
	}
	
	Scalar PetscBackend::reduce(const PETScVector &v, const Min &)
	{
		Scalar x = 0.0;
		VecMin(v.implementation(), nullptr, &x);
		return x;
	}

	template<class Operation>
	inline static Scalar generic_local_reduce(const PETScMatrix &m, const Scalar &init_value, const Operation &op)
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

	Scalar PetscBackend::reduce(const PETScMatrix &m, const Min &op)
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
			VecCreateMPI(comm, lrows, grows, &v);
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

	Scalar PetscBackend::reduce(const PETScVector &v, const Max &)
	{
		Scalar x = -std::numeric_limits<Scalar>::max();
		VecMax(v.implementation(), nullptr, &x);
		return x;
	}

	Scalar PetscBackend::reduce(const PETScMatrix &m, const Max &op)
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
			VecCreateMPI(comm, lrows, grows, &v);
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
	void PetscBackend::diag(PETScVector &vec, const PETScMatrix &mat)
	{
		using std::max;
		PetscInt globalRows, globalColumns;
		
		PetscErrorCode err = MatGetSize(mat.implementation(), &globalRows, &globalColumns);
		PetscInt localRows, localColumns;
		err = MatGetLocalSize(mat.implementation(), &localRows, &localColumns);
		
		PetscInt lenGlobal = max(globalRows, globalColumns);
		PetscInt lenLocal  = max(localRows, localColumns);
		resize(vec, {lenLocal}, {lenGlobal});
		err = MatGetDiagonal(mat.implementation(), vec.implementation());
		PETScError::Check(err);
	}
	
	void PetscBackend::diag(PETScSparseMatrix &mat, const PETScVector &vec)
	{
		// I do not think, this is needed
		// because, doesnt run in parallel properly
		// u can not change size of matrix after it was initialized...
		// or at least it seems like ...
		//!!! FIXME needs to check if matrix has been allocated already
		
		MPI_Comm comm = mat.communicator();
		MatDestroy(&mat.implementation());
		MatCreate(comm, &mat.implementation());

		PetscInt local_size, global_size;
		VecGetLocalSize(vec.implementation(), &local_size);
		VecGetSize(vec.implementation(), &global_size);
		
		if(mpi_world_size() == 1) {
			MatSetType(mat.implementation(), MATSEQAIJ);
			//
			// MatSEQAIJSetPreallocation(mat.implementation(), 1, NULL, NULL, NULL);
			MatSetSizes(mat.implementation(), local_size, local_size, global_size, global_size);
			MatSetUp(mat.implementation());
		} else {
			MatSetType(mat.implementation(), MATMPIAIJ);
			MatSetSizes(mat.implementation(), local_size, local_size, global_size, global_size);
			MatMPIAIJSetPreallocation(mat.implementation(), 1, NULL, 0, NULL);

			MatSetOption(mat.implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE);
			MatSetOption(mat.implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
			MatSetOption(mat.implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);
		}
		
		
		PetscInt err = MatDiagonalSet( mat.implementation(), vec.implementation(), INSERT_VALUES);
		PETScError::Check(err);
		
	}

	void PetscBackend::diag(PETScMatrix &out, const PETScMatrix &in)
	{
		PETScVector vec;
		diag(vec, in);
		diag(out, vec);
	}

	void PetscBackend::diag(PETScMatrix &mat, const PETScVector &vec)
	{
		// I do not think, this is needed
		// because, doesnt run in parallel properly
		// u can not change size of matrix after it was initialized...
		// or at least it seems like ...
		//!!! FIXME needs to check if matrix has been allocated already
		
		MPI_Comm comm = mat.communicator();
		MatDestroy(&mat.implementation());
		MatCreate(comm, &mat.implementation());
		
		PetscInt local_size, global_size;
		VecGetLocalSize(vec.implementation(), &local_size);
		VecGetSize(vec.implementation(), &global_size);
		
		if(mpi_world_size() == 1) {
			MatSetType(mat.implementation(), MATSEQDENSE);
			//
			// MatSEQAIJSetPreallocation(mat.implementation(), 1, NULL, NULL, NULL);
			MatSetSizes(mat.implementation(), local_size, local_size, global_size, global_size);
			MatSetUp(mat.implementation());
		} else {
			MatSetType(mat.implementation(), MATMPIDENSE);
			MatSetSizes(mat.implementation(), local_size, local_size, global_size, global_size);
			// MatMPIAIJSetPreallocation(mat.implementation(), 1, NULL, 0, NULL);
			MatSetUp(mat.implementation());
		}
		
		
		PetscInt err = MatDiagonalSet( mat.implementation(), vec.implementation(), INSERT_VALUES);
		PETScError::Check(err);
		
	}

	void PetscBackend::mat_diag_shift(PETScMatrix &left, const Scalar diag_factor)
	{
		PETScError::Check(MatShift(left.implementation(), diag_factor));
	}
	
	bool PetscBackend::compare(const Vector &left, const Vector &right, const ApproxEqual &comp) {
		PETScVector diff;
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
		MatCreateDense(comm, l_range.extent(), PETSC_DECIDE, result_size.get(0), result_size.get(1), NULL,
					   &result.implementation());
		
		
		for(SizeType i = l_range.begin(); i != l_range.end(); ++i) {
			const Scalar l_value = left_array[i-l_range.begin()];
			
			for(SizeType j = 0; j != n; ++j) {
				const Scalar r_value = recvbuf.at(j);
				
				MatSetValue(result.implementation(), i, j, l_value*r_value, INSERT_VALUES);
			}
		}
		
		MatAssemblyBegin(result.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(result.implementation(), MAT_FINAL_ASSEMBLY);
		VecRestoreArrayRead(left.implementation(), &left_array);
	}
	

	void PetscBackend::assign_to_range(PETScMatrix & /*left*/, const PETScMatrix &/*right*/, const Range &/*global_row_range*/,
									 const Range &/*global_col_range*/)
	{
		
		assert(false); //TODO
		
	}
	
	void PetscBackend::assign_to_range( PETScMatrix &/*left*/, const Identity &/**/, const Range &/*global_row_range*/,
									 const Range &/*global_col_range*/) {
		assert(false); //TODO
	}
	
	void PetscBackend::assign_to_range( PETScVector &/*left*/, const PETScVector &/*right*/, const Range &/*global_row_range*/,
									 const Range &/*global_col_range*/) {
		assert(false); //TODO
	}
	
	void PetscBackend::vec_to_mat(Matrix &m, const Vector & v, const bool transpose) {
		
		PetscInt n;
		VecGetSize(v.implementation(), &n);
		Range r = range(v);
		
		Matrix mtr;
		
		Matrix &work = transpose ? mtr : m;
		build(work, Size({n, 1}), Values<Scalar>(0));
		
		
		
		Scalar zero = 0;
		for (PetscInt i = r.begin(); i < r.end(); ++i) {
			Scalar val = get(v, i);
			set(work, i, zero, val);
		}
		
		MatAssemblyBegin(work.implementation(), MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(work.implementation(), MAT_FINAL_ASSEMBLY);
		
		if (transpose) {
			MatTranspose(mtr.implementation(), MAT_INITIAL_MATRIX, &m.implementation());
		}

	}
	
	// this is not very correct yet ...
	void PetscBackend::build_local_diag_block(PETScSerialSparseMatrix &left, const PETScSparseMatrix &right)
	{
		Mat M;
		MatGetDiagonalBlock(right.implementation(), &M);
		convert(M, left);
	}
	
	void PetscBackend::triple_product_ptap(PETScMatrix & result, const PETScMatrix & A, const PETScMatrix &P)
	{
		if(&result.implementation() != &A.implementation() && &result.implementation() != &P.implementation()) {
			MatDestroy(&result.implementation());
		} else {
			std::cerr << "[Warning] not handled case in triple_product_ptap" << std::endl;
		}

		MatPtAP(A.implementation(), P.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()); 
	}
	


	void PetscBackend::triple_product(PETScMatrix & result, const PETScMatrix &A, const PETScMatrix &B, const PETScMatrix &C)
	{
		if(&result.implementation() != &A.implementation() && &result.implementation() != &B.implementation() && &result.implementation() != &C.implementation()) {
			MatDestroy(&result.implementation());
		}
		
		MatMatMatMult(A.implementation(), B.implementation(), C.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()); 
	}

	bool PetscBackend::is_nan_or_inf(const PETScVector &X)
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

	bool PetscBackend::is_nan_or_inf(const PETScMatrix &m)
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

		MPI_Comm comm = PetscObjectComm((PetscObject) m.implementation());
		MPI_Allreduce(MPI_IN_PLACE, &has_nan, 1, MPI_INT, MPI_MAX, comm);
		return has_nan > 0;
	}

	// redistribution of local sizes of vector
	// TODO: can be done also for matrices 
	//       can be done also based on local sizes, no provided vector
	void PetscBackend::build_local_redistribute(PETScVector &result, const PETScVector &x_from, const PETScVector &shape_vec)
	{
		if(mpi_world_size() == 1) 
        {	
        	result = x_from; 
        } 
        else
        {               
            PETScVector x_to = shape_vec; 
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


	void PetscBackend::axpy(PETScVector &y, const Scalar alpha, const PETScVector &x) 
	{
		VecAXPY(y.implementation(), alpha, x.implementation());
	}

	void PetscBackend::axpy(PETScMatrix &y, const Scalar alpha, const PETScMatrix &x)
	{
		MatAXPY(y.implementation(), alpha, x.implementation(), DIFFERENT_NONZERO_PATTERN);
	}

	void PetscBackend::axpby(Vector &y, const Scalar alpha, const Vector &x, const Scalar &beta)
	{
		VecAXPBY(y.implementation(), alpha, beta, x.implementation());
	}
	


	void PetscBackend::diag_scale_right(Matrix &result, const Matrix &m, const Vector &diag)
	{
		assign(result, m);
		PETScError::Check( MatDiagonalScale(result.implementation(), nullptr, diag.implementation()) );
	}

	void PetscBackend::diag_scale_left(Matrix &result, const Vector &diag, const Matrix &m)
	{
		assign(result, m);
		PETScError::Check( MatDiagonalScale(result.implementation(), diag.implementation(), nullptr) );
	}

	void PetscBackend::diag_scale_left(Vector &result, const Vector &diag, const Vector &m)
	{
        Size gs, ls;
        size(m, gs);
        local_size(m, ls);
        resize(result, ls, gs);
        PETScError::Check( VecPointwiseMult(diag.implementation(), m.implementation(), result.implementation()) );
	}

	Scalar PetscBackend::trace(const Matrix &mat)
	{
		Scalar ret;
		MatGetTrace(mat.implementation(), &ret);
		return ret;
	}

	template<class Operation>
	inline static void generic_col_reduce(
		PETScVector &result,
		const PETScMatrix &mat, 
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
		VecCreateMPI(mat.communicator(), local_r, global_r, &result.implementation());
		

		// auto fun = Operation::template Fun<Scalar>();

		MatGetOwnershipRange(mat.implementation(), &r_begin, &r_end);
		
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

		VecAssemblyBegin(result.implementation());
		VecAssemblyEnd(result.implementation());
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Plus &, const int dim)
	{
		PETScVector rowSum; //FIXME initialize
		Size gs, ls;
		size(mat, gs);
		local_size(mat, ls);
		resize(result, ls, gs);

		if(dim == 1) {
			MatGetRowSum(mat.implementation(), result.implementation());
		} else {
			//FIXME implement own
			assert(false && "not available in pestsc");
		}
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Min &op, const int dim)
	{
		// Scalar x;
		PetscInt grows, gcols;
		MatGetSize(mat.implementation(), &grows, &gcols);

		int comm_size = 0;
		MPI_Comm_size(mat.communicator(), &comm_size);
		bool is_serial = comm_size == 1;

		if (dim == 1) {
			if(is_serial) {
				VecDestroy(&result.implementation());
				VecCreateMPI(mat.communicator(), PETSC_DECIDE, grows, &result.implementation());
				MatGetRowMin(mat.implementation(), result.implementation(), nullptr);
			} else {
				generic_col_reduce(
					result,
					mat, 
					std::numeric_limits<Scalar>::max(), 
					op
					);
			}
		} else {
			//FIXME implement own
			// VecDestroy(&result.implementation());
			// VecCreateMPI(mat.communicator(), PETSC_DECIDE, gcols, &result.implementation());
			assert(false && "not available in petsc");
		}
	}

	void PetscBackend::apply_tensor_reduce(Vector &result, const Matrix &mat, const Max &op, const int dim)
	{
		// Scalar x;
		PetscInt grows, gcols;
		MatGetSize(mat.implementation(), &grows, &gcols);

		int comm_size = 0;
		MPI_Comm_size(mat.communicator(), &comm_size);
		bool is_serial = comm_size == 1;

		if (dim == 1) {
			if(is_serial) {
			VecDestroy(&result.implementation());
			VecCreateMPI(mat.communicator(), PETSC_DECIDE, grows, &result.implementation());
			MatGetRowMax(mat.implementation(), result.implementation(), nullptr);
			} else {
				generic_col_reduce(
					result,
					mat, 
					-std::numeric_limits<Scalar>::max(), 
					op);
			}
		} else {
			//FIXME implement own
			// VecDestroy(&result.implementation());
			// VecCreateMPI(mat.communicator(), PETSC_DECIDE, gcols, &result.implementation());
			assert(false && "not available in petsc");
		}

	}

	void PetscBackend::inverse(Matrix &result, const Matrix &mat)
	{
		PETScMatrix I;

		Size s;
		size(mat, s);

		build(I, s, Identity());
		build(result, s, Zeros());

		IS isr, isc;
		MatFactorInfo info;

		Matrix L;
		assign(L, mat);

		MatGetOrdering(L.implementation(), MATORDERINGNATURAL, &isr, &isc);
		
		MatLUFactor(L.implementation(), isr, isc, &info); 
		PETScError::Check( MatMatSolve(L.implementation(), I.implementation(), result.implementation()) );
		ISDestroy(&isr);
		ISDestroy(&isc);
	}
}

