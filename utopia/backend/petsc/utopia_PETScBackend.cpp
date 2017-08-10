#include "utopia_PETScBackend.hpp"

#include "petscmat.h"
#include "petscvec.h"

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


		inline const Mat &left() const
		{
			return left_;
		}

		inline const Mat &right() const
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
			// std::cout << type_str.substr(start, 5) << std::endl;
			return !(type_str.substr(start, 5) == "dense");
		}

		Mat left_;
		Mat right_;
		bool must_destroy_left_;
		bool must_destroy_right_;
	};

	void PETScBackend::resize(const Size &s, PETScMatrix &mat)
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

	void PETScBackend::resize(const Size &local_s, const Size &global_s, PETScMatrix &mat)
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

	void PETScBackend::resize(const Size &s, PETScVector &vec)
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


	void PETScBackend::resize(const Size &s_local, const Size &s_global, PETScVector &vec)
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

	bool PETScBackend::transpose(const PETScMatrix &mat, PETScMatrix &result)
	{
		if(&mat != &result) {
			MatDestroy(&result.implementation());
		} else {
			assert(false);
		}

		return PETScError::Check(MatTranspose(mat.implementation(), MAT_INITIAL_MATRIX, &result.implementation()));
	}
	
	void PETScBackend::clear(PETScMatrix &mat) const
	{
		MatDestroy(&mat.implementation());
		MatCreate(mat.communicator(), &mat.implementation());
	}
	
	
	bool PETScBackend::convert(Vec vec, PETScVector &wrapper)
	{
		MPI_Comm comm = PetscObjectComm((PetscObject) vec);

		PetscInt n, local_n;
		VecGetSize(vec, &n);
		VecGetLocalSize(vec, &local_n);
		VecDestroy(&wrapper.implementation());
		VecCreateMPI(comm, local_n, n, &wrapper.implementation());
		VecCopy(vec, wrapper.implementation());
		return true;
	}
	
	bool PETScBackend::convert(Mat mat, PETScMatrix &wrapper)
	{
		MPI_Comm comm = PetscObjectComm((PetscObject) mat);

		PetscInt r, c, local_r, local_c;
		MatGetSize(mat, &r, &c);
		MatGetLocalSize(mat, &local_r, &local_c);
		
		MatDestroy(&wrapper.implementation());
		MatCreateDense(comm, local_r, local_c, r, c, NULL,
					   &wrapper.implementation());
		MatCopy(mat, wrapper.implementation(), SAME_NONZERO_PATTERN);
		return true;
	}
	
	bool PETScBackend::convert(Mat mat, PETScSparseMatrix &wrapper)
	{
		//            MatCreateAIJ(m.communicator(), PETSC_DECIDE, PETSC_DECIDE, rows, cols, 1, PETSC_NULL,
		//                         1 /*Only because otherwise petsc crashes*/, PETSC_NULL, &m.implementation());
		MatDestroy(&wrapper.implementation());
		MatDuplicate(mat, MAT_COPY_VALUES, &wrapper.implementation());
		return true;
	}
	
	bool PETScBackend::convert(PETScVector &wrapper, Vec vec)
	{
		VecCopy(wrapper.implementation(), vec);
		return true;
	}
	
	bool PETScBackend::convert(PETScMatrix &wrapper, Mat mat)
	{
		MatCopy(wrapper.implementation(), mat, SAME_NONZERO_PATTERN);
		return true;
	}
	
	bool PETScBackend::convert(PETScSparseMatrix &wrapper, Mat mat)
	{
		MatCopy(wrapper.implementation(), mat, SAME_NONZERO_PATTERN);
		return true;
	}

	bool PETScBackend::wrap(Mat mat, PETScSparseMatrix &wrapper)
	{
		wrapper.wrap(mat);
		return true;
	}


	bool PETScBackend::wrap(Vec vec, PETScVector &wrapper)
	{
		// wrapper.wrap(vec);
		return true;
	}
	
	
	
	Range PETScBackend::range(const PETScVector &v) 
	{
		PetscInt rbegin, rend;
		VecGetOwnershipRange(v.implementation(), &rbegin, &rend);
		return Range(rbegin, rend);
	}
	
	Range PETScBackend::rowRange(const PETScMatrix &m) 
	{
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		assert(Range(rbegin, rend).valid());
		return Range(rbegin, rend);
	}
	
	Range PETScBackend::colRange(const PETScMatrix &m) 
	{
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		return Range(0, gcols);
	}
	
	bool PETScBackend::size(const PETScMatrix &m, Size &size) const 
	{
		PetscInt grows, gcols;
		size.setDims(2);
		MatGetSize(m.implementation(), &grows, &gcols);
		size.set(0, grows);
		size.set(1, gcols);
		return true;
	}
	
	bool PETScBackend::size(const PETScVector &m, Size &size) const 
	{
		PetscInt n;
		size.setDims(1);
		VecGetSize(m.implementation(), &n);
		size.set(0, n);
		return true;
	}
	
	bool PETScBackend::local_size(const PETScVector &m, Size &size) const 
	{
		PetscInt n;
		size.setDims(1);
		VecGetLocalSize(m.implementation(), &n);
		size.set(0, n);
		return true;
	}

	bool PETScBackend::local_size(const PETScMatrix &mat, Size &size) const
	{
		PetscInt n, m;
		size.setDims(2);
		MatGetLocalSize(mat.implementation(), &n, &m);
		size.set(0, n);
		size.set(1, m);
		return true;
	}
	
	void PETScBackend::assignFromRange(PETScMatrix &left, const PETScMatrix &right, const Range &globalRowRange,
									   const Range &globalColRange) {
		assert(!globalRowRange.empty());
		
		Mat &li = left.implementation();
		//            const Mat &ri = left.implementation();

		MPI_Comm comm = right.communicator();
		MatDestroy(&left.implementation());
		
		MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, globalRowRange.extent(),
					   globalColRange.extent(), PETSC_NULL, &li);
		
		
		Range rr = rowRange(right).intersect(globalRowRange);
		Range cr = colRange(right).intersect(globalColRange);
		
		// rr might be invalid !!! so use < to compare with range
		//            if(!rr.valid()) std::cout << rr << std::endl;
		
		
		//use MAT_FLUSH_ASSEMBLY to not make petsc crash if nothing is added
		MatAssemblyBegin(li, MAT_FLUSH_ASSEMBLY);
		
		PetscInt r = 0, c = 0;
		for (PetscInt rIt = rr.begin(); rIt < rr.end(); ++rIt) {
			r = rIt - globalRowRange.begin();
			for (PetscInt cIt = cr.begin(); cIt < cr.end(); ++cIt) {
				c = cIt - globalColRange.begin();
				set(left, r, c, get(right, rIt, cIt));
			}
		}
		
		MatAssemblyEnd(li, MAT_FLUSH_ASSEMBLY);
		
		//To finalize the matrix
		MatAssemblyBegin(li, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(li, MAT_FINAL_ASSEMBLY);
	}
	
	
	// read matrix
	bool PETScBackend::read(const std::string &path, PETScMatrix &Mat_A)
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
	bool PETScBackend::write(const std::string &path, const PETScMatrix &Mat_A)
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
	bool PETScBackend::write(const std::string &path, const PETScVector &Vec_A)
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
	bool PETScBackend::monitor(const long & it, PETScMatrix &Mat_A)
	{
		PetscViewer viewer_hessian = nullptr;
		if(it==0)
		{
          	// log stifness
          	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "log_hessian.m" ,&viewer_hessian);  
          	PetscViewerPushFormat(viewer_hessian,PETSC_VIEWER_ASCII_MATLAB); 
        }   
    	Mat A; 
    	MatDuplicate(Mat_A.implementation(), MAT_COPY_VALUES,  &A); 
    	MatCopy(Mat_A.implementation(), A, SAME_NONZERO_PATTERN); 

	    PetscObjectSetName((PetscObject)A, "H");
	    MatView(A, viewer_hessian); 
		//PetscViewerDestroy(&viewer_hessian);
    
		return true; 
	}
	
	// monitor for cyrill  - maybe u need to adjust something here ...
	bool PETScBackend::monitor(const long & it, PETScVector &Vec_A)
	{
		PetscViewer viewer_iterates = nullptr;
		if(it==0)
		{
          	// log iterates
          	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "log_iterate.m", &viewer_iterates);  
          	PetscViewerPushFormat(viewer_iterates,PETSC_VIEWER_ASCII_MATLAB); 
        }
        Vec iterates; 	        
        VecDuplicate(Vec_A.implementation(), &iterates);
  		VecCopy(Vec_A.implementation(),iterates);
  		PetscObjectSetName((PetscObject)iterates,"it");
  		VecView(iterates, viewer_iterates); 

      	//PetscViewerDestroy(&viewer_iterates);
		return true; 
	}

	PetscScalar PETScBackend::get_global_nnz(PETScMatrix &Mat_A)
	{	
		MatInfo        info;
		MatGetInfo(Mat_A.implementation(), MAT_GLOBAL_SUM, &info); 
		return info.nz_used; 
	}


	PetscScalar PETScBackend::get_local_nnz(PETScMatrix &Mat_A)
	{	
		MatInfo        info;
		MatGetInfo(Mat_A.implementation(), MAT_LOCAL, &info); 
		return info.nz_used; 
	}



	bool PETScBackend::set_zero_rows(PETScMatrix &Mat_A, const std::vector<int> &index)
	{
		return PETScError::Check(MatZeroRows(Mat_A.implementation(), index.size(), &index[0], 1.0, NULL, NULL));
	}

	bool  PETScBackend::apply_BC_to_system(PETScMatrix & A, PETScVector& x, PETScVector& rhs, const std::vector<int> &index)
	{
		return PETScError::Check(MatZeroRows(A.implementation(), index.size(), &index[0], 1.0, x.implementation(), rhs.implementation())); 
	}

	// read vector
	bool PETScBackend::read(const std::string &path, PETScVector &Vec_A)
	{
		MPI_Comm comm = Vec_A.communicator();

		Vec &A = Vec_A.implementation();
		VecDestroy(&A);

		PetscViewer fd;
		PetscViewerBinaryOpen(comm, path.c_str(), FILE_MODE_READ, &fd);
		VecCreate(comm,&A);
		bool status;
		status =  PETScError::Check( VecLoad(A,fd));
		PetscViewerDestroy(&fd);
		return status;
	}
	
	void PETScBackend::assignFromRange(PETScVector &left, const PETScVector &right, const Range &globalRowRange,
									   const Range & /*globalColRange */) {
		assert(!globalRowRange.empty());
		
		Vec &li = left.implementation();
		//const Vec &ri = left.implementation();
		
		// get range of vector
		Range rr = range(right).intersect(globalRowRange);
		
		// create vector
		//VecCreate(left.communicator(), &li);
		// TODO Check VecCreate on already exisiting vec
		VecSetType(li, VECMPI);
		
		// set vector size
		VecSetSizes(li, PETSC_DECIDE, rr.extent());
		//left.initialize(PETSC_DECIDE, rr.extent());
		
		VecAssemblyBegin(li);
		
		PetscInt r = 0;
		for (PetscInt rIt = rr.begin(); rIt < rr.end(); ++rIt) {
			r = rIt - globalRowRange.begin();
			set(left, r, get(right, rIt));
		}
		
		//VecAssemblyEnd(li);
		
		// Finialize
		//VecAssemblyBegin(li);
		VecAssemblyEnd(li);
		
		//assert(false); //TODO
	}
	
	bool PETScBackend::scal(const PetscScalar scaleFactor, const PETScMatrix &left , PETScMatrix &right) {
		if(&left == &right) {
			MatScale(right.implementation(), scaleFactor);
		} else {
			assert(false); // FIXME
			return false;
		}

		return true;
	}
	
	
	bool PETScBackend::gemm(const PetscScalar /*alpha*/, const PETScMatrix &/*left*/, const PETScMatrix &/*right*/, const PetscScalar /* beta*/,
							PETScMatrix & /*result*/) {
		assert(false); //TODO
		return false;
	}

	bool PETScBackend::gemm(const PetscScalar alpha, const PETScMatrix &left, const PETScMatrix &right,
							bool transpose_left, bool transpose_right, const PetscScalar beta, PETScMatrix &result) {
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
				ok = transpose(temp, result);
			}
		}

		return ok;
	}

	bool PETScBackend::gemm(const double scaleFactor, const PETScMatrix &left, const Vector &right, bool transpose_left,  bool transpose_right, const double beta, Vector &result)
	{
		 assert(!transpose_right);
		
		 Size gs, ls;
		 size(left, gs);
		 local_size(left, ls);

		 VecSetType(result.implementation(), VECMPI);
		 
		 if(transpose_left) {
		 	VecSetSizes(result.implementation(), ls.get(1), gs.get(1));
		 	return PETScError::Check( MatMultTranspose(left.implementation(), right.implementation(), result.implementation()) );
		 } else {
		 	VecSetSizes(result.implementation(), ls.get(0), gs.get(0));
		 	return PETScError::Check( MatMult(left.implementation(), right.implementation(), result.implementation()) );
		 }
	}
	
	bool PETScBackend::gemm(const PetscScalar /*alpha */, const PETScVector &/*left*/, const PETScMatrix &/*right*/,
							bool /*transpose_left*/, bool /*transpose_right*/, const PetscScalar /*beta*/, PETScMatrix &/*result*/) {
		assert(false); //TODO
		return false;
	}
	
	bool PETScBackend::scal(const PetscScalar /*scale_factor*/, const PETScVector &/*v*/, PETScVector &/*result */) {
		assert(false); //TODO
		return false;
	}
	
	void PETScBackend::build(PETScMatrix &m, const Size &size, const Identity &) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		//FIXME use this: MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag,Vec x,Vec b)
		MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), NULL,
					   &m.implementation());
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		
		rend = PetscMin(rend, gcols);
		for (PetscInt i = rbegin; i < rend; ++i) {
			MatSetValue(m.implementation(), i, i, 1, INSERT_VALUES);
		}
		
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::build(PETScSparseMatrix &m, const Size &size, const Identity &) {
		
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);
		//PetscInt n = std::min(rows, cols);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, rows, cols, 1, PETSC_NULL,
					 1 /*Only because otherwise petsc crashes*/, PETSC_NULL, &m.implementation());
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		rend = PetscMin(rend, gcols);
		
		for (PetscInt i = rbegin; i < rend; ++i) {
			const PetscInt offset = i;
			const PetscScalar val = 1;
			MatSetValues(m.implementation(), 1, &offset, 1, &offset, &val, INSERT_VALUES);
		}
		
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::build(PETScMatrix &m, const Size &size, const LocalIdentity &) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		//FIXME use this: MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag,Vec x,Vec b)
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
		
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		
		PetscInt extent = PetscMin(rend-rbegin, PetscMin(grows, gcols));
		
		
		for (PetscInt i = 0; i < extent; ++i) {
			MatSetValue(m.implementation(), rbegin+i, receive_buffer + i, 1, INSERT_VALUES);
		}
		
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::build(PETScSparseMatrix &m, const Size &size, const LocalIdentity &) {
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);
		//PetscInt n = std::min(rows, cols);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, rows, cols, PETSC_DETERMINE, PETSC_DETERMINE, 1, PETSC_NULL,
					 1 /*Only because otherwise petsc crashes*/, PETSC_NULL, &m.implementation());
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		//PetscInt grows, gcols;
		//MatGetSize(m.implementation(), &grows, &gcols);
		
		unsigned long send_buffer = cols;
		unsigned long receive_buffer = 0;
		
		MPI_Exscan(&send_buffer, &receive_buffer, 1, MPI_UNSIGNED_LONG ,
				   MPI_SUM, m.communicator());
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt extent = PetscMin(rend-rbegin, PetscMin(rows, cols));
		
		for (PetscInt i = 0; i < extent; ++i) {
			const PetscInt row_offset = rbegin+i;
			const PetscInt col_offset = receive_buffer+i;
			const PetscScalar val = 1;
			MatSetValues(m.implementation(), 1, &row_offset, 1, &col_offset, &val, INSERT_VALUES);
		}
		
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		
	}
	
	void PETScBackend::build(PETScSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz) {
		
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, rows, cols,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to local entries*/, PETSC_NULL,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to remote entries*/, PETSC_NULL,
					 &m.implementation());
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		
		
	}
	
	void PETScBackend::build(PETScSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz) {
		
		//Sparse id
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, rows, cols, PETSC_DETERMINE, PETSC_DETERMINE,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to local entries*/, PETSC_NULL,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to remote entries*/, PETSC_NULL,
					 &m.implementation());
		
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		
	}
	
	void PETScBackend::build(PETScSparseMatrix &m, const Size &size, const LocalRowNNZ<PetscInt> &nnz)
	{
		PetscInt rows = size.get(0);
		PetscInt cols = size.get(1);

		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());
		
		MatCreateAIJ(comm, rows, PETSC_DECIDE, PETSC_DETERMINE, cols,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to local entries*/, PETSC_NULL,
					 PetscMax(nnz.nnz(), 1) /*n DOF connected to remote entries*/, PETSC_NULL,
					 &m.implementation());
		
		
		MatSetOption(m.implementation(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	}
	
	/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
	void PETScBackend::build(PETScMatrix  &m, const Size &size, const LocalNNZ<PetscInt> & /*nnz */)
	{
		build(m, size, LocalValues<PetscScalar>(0));
	}
	
	/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
	void PETScBackend::build(PETScMatrix  &m, const Size &size, const NNZ<PetscInt> &/*nnz*/)
	{
		build(m, size, Zeros());
	}
	
	void PETScBackend::build(PETScMatrix &m, const Size &size, const Zeros &)
	{
		build(m, size, Values<PetscScalar>(0));
	}
	
	void PETScBackend::build(PETScVector &v, const Size &size, const Zeros &)
	{
		build(v, size, Values<PetscScalar>(0));
	}
	
	void PETScBackend::build(PETScMatrix &m, const Size &size, const LocalZeros &)
	{
		build(m, size, LocalValues<PetscScalar>(0));
	}
	
	void PETScBackend::build(PETScVector &v, const Size &size, const LocalZeros &)
	{
		build(v, size, LocalValues<PetscScalar>(0));
	}
	
	void PETScBackend::build(PETScMatrix &m, const Size &size, const Values<PetscScalar> &values) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, size.get(0), size.get(1), NULL,
					   &m.implementation());
		
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		const PetscScalar v = values.value();
		for (PetscInt i = rbegin; i < rend; ++i) {
			for (PetscInt j = 0; j < gcols; ++j) {
				MatSetValue(m.implementation(), i, j, v, INSERT_VALUES);
			}
		}
		
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::build(PETScVector &v, const Size &local_size, const Size &&global_size, const Values<PetscScalar> &values)
	{
		Size v_local_size;
		if(v.isInitialized()) {
			this->local_size(v, v_local_size);
		}	

		if(!v.isInitialized() || v_local_size.get(0) != local_size.get(0)) {
			VecDestroy(&v.implementation());
			VecCreateMPI(v.communicator(), local_size.get(0), global_size.get(0), &v.implementation());
		}
 
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());
	}
	
	void PETScBackend::build(PETScVector &v, const Size &size, const Values<PetscScalar> &values) {

		VecDestroy(&v.implementation());
		VecCreateMPI(v.communicator(), PETSC_DECIDE, size.get(0), &v.implementation());
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());
	}
	
	void PETScBackend::build(PETScMatrix &m, const Size &size, const LocalValues<PetscScalar> &values) {
		MPI_Comm comm = m.communicator();
		MatDestroy(&m.implementation());

		MatCreateDense(comm, size.get(0), size.get(1), PETSC_DETERMINE, PETSC_DETERMINE, NULL,
					   &m.implementation());
		
		//TODO check if it can be simplified using known information.
		
		PetscInt rbegin, rend;
		MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
		
		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		
		
		MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);
		
		
		const PetscScalar v = values.value();
		for (PetscInt i = rbegin; i < rend; ++i) {
			for (PetscInt j = 0; j < gcols; ++j) {
				MatSetValue(m.implementation(), i, j, v, INSERT_VALUES);
			}
		}
		
		MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::build(PETScVector &v, const Size &size, const LocalValues<PetscScalar> &values) {
		MPI_Comm comm = v.communicator();
		VecDestroy(&v.implementation());

		VecCreateMPI(comm, size.get(0), PETSC_DETERMINE, &v.implementation());
		
		VecSet(v.implementation(), values.value());
		VecAssemblyBegin(v.implementation());
		VecAssemblyEnd(v.implementation());
	}
	
	void PETScBackend::add(PETScVector &v, const PetscInt index, PetscScalar value)
	{
		VecSetValues(v.implementation(), 1, &index, &value, ADD_VALUES);
	}
	
	void PETScBackend::add(PETScMatrix &m, const PetscInt row, const PetscInt col, PetscScalar value)
	{
		MatSetValues(m.implementation(), 1, &row, 1, &col, &value, ADD_VALUES);
	}
	
	void PETScBackend::set(PETScVector &v, const PetscInt index, PetscScalar value) {
		VecSetValues(v.implementation(), 1, &index, &value, INSERT_VALUES);
	}
	
	void PETScBackend::set(PETScVector &v, const std::vector<PetscInt> indices, const std::vector<PetscScalar> values) {
		assert(indices.size() == values.size());
		VecSetValues(v.implementation(), indices.size(), &indices[0], &values[0], INSERT_VALUES);
	}
	
	
	void PETScBackend::set(PETScMatrix &v, const PetscInt row, const PetscInt col, PetscScalar value) {
		MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
	}
	
	void PETScBackend::set(PETScSparseMatrix &v, const PetscInt row, const PetscInt col, PetscScalar value) {
		MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
	}
	
	void PETScBackend::writeLock(PETScVector &vec) {
		vec.assemblyBegin();
	}
	
	void PETScBackend::writeUnlock(PETScVector &vec) {
		vec.assemblyEnd();
	}
	
	void PETScBackend::writeLock(const PETScMatrix &mat) {
		MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::writeUnlock(const PETScMatrix &mat) {
		MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::writeLock(const PETScSparseMatrix &mat) {
		MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::writeUnlock(const PETScSparseMatrix &mat) {
		MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
	}
	
	void PETScBackend::set(PETScMatrix &v, const std::vector<PetscInt> rows, const std::vector<PetscInt> cols,
						   const std::vector<PetscScalar> values) {
		assert(rows.size() == values.size());
		assert(cols.size() == values.size());
		//FIXME use this instead of the loop below (PETSC has buggy behavior though)
		//            MatSetValues(v.implementation(), static_cast<PetscInt>(rows.size()), &rows[0],
		//                                             static_cast<PetscInt>(cols.size()), &cols[0],
		//                                             &values[0], INSERT_VALUES);
		
		for(std::vector<PetscInt>::size_type i = 0; i != rows.size(); ++i) {
			set(v, rows[i], cols[i], values[i]);
		}
	}
	
	PetscScalar PETScBackend::get(const PETScVector &v, const PetscInt index) {
		PetscScalar value;
		VecGetValues(v.implementation(), 1, &index, &value);
		return value;
	}
	
	PetscScalar PETScBackend::get(const PETScMatrix &v, const PetscInt row, const PetscInt col) {
		PetscScalar value;
		MatGetValues(v.implementation(), 1, &row, 1, &col, &value);
		return value;
	}
	
	bool PETScBackend::apply(const PETScMatrix &left, const PETScVector &right, const Multiplies &, PETScVector &result) {
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
		
		return true;
	}
	
	bool PETScBackend::apply(const PETScMatrix &left, const PETScMatrix &right, const Multiplies &, PETScMatrix &result)
	{

		if(&right.implementation() != &result.implementation() || &left.implementation() !=  &result.implementation()) {
			MatDestroy(&result.implementation());
		} else {
			assert(false);
		}
 

		MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
		return true;
	}

	bool PETScBackend::mat_mult_add(const PETScMatrix &m, const PETScVector &right, const PETScVector &left, PETScVector &result)
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
		return true;
	}

	bool PETScBackend::mat_multT_add(const PETScMatrix &m, const PETScVector &right, const PETScVector &left, PETScVector &result)
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
		return true;
	}

	bool PETScBackend::apply(const Vector &vec, const Abs &, Vector &result)
	{
		if(&vec != &result) {
			// bool keep_size = false;
			Size l_s, g_s;

			local_size(vec, l_s);
			size(vec, g_s);

			resize(l_s, g_s, result);
			PETScError::Check( VecCopy(vec.implementation(), result.implementation()) );
			
		}	

		VecAbs(result.implementation());
		return true;
	}

	// bool PETScBackend::apply(const Matrix &mat, const Abs &, Matrix &result)
	// {
	// 	if(&mat != &result) {
	// 		bool keep_size = false;
	// 		Size l_s, g_s;

	// 		local_size(mat, l_s);
	// 		size(mat, g_s);

	// 		resize(l_s, g_s, result);
	// 		PETScError::Check( MatCopy(mat.implementation(), result.implementation(), SAME_NONZERO_PATTERN) );
	// 	}	

	// 	MatAbs(result.implementation());
	// 	return true;
	// }


	static void allocate_apply_vec(const PETScVector &left, const PETScVector &right, PETScVector &result)
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
	
	bool PETScBackend::apply(const PETScVector &left, const PETScVector &right, const EMultiplies &, PETScVector &result) {
		allocate_apply_vec(left, right, result);
		return PETScError::Check(VecPointwiseMult(result.implementation(), left.implementation(), right.implementation()));
	}
	
	bool PETScBackend::apply(const PETScVector &left, const PETScVector &right, const Divides &, PETScVector &result) {
		allocate_apply_vec(left, right, result);
		return PETScError::Check(VecPointwiseDivide(result.implementation(), left.implementation(), right.implementation()));
	}
	
	// reciprocal
	bool PETScBackend::apply(const PETScVector &vec, const Reciprocal<double> &reciprocal, PETScVector &result)
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

		return PETScError::Check(err);
	}
	
	PetscScalar PETScBackend::norm2(const PETScVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_2, &val);
		return val;
	}
	
	PetscScalar PETScBackend::norm2(const Matrix &m) {
		PetscReal val;
		MatNorm(m.implementation(), NORM_FROBENIUS, &val);
		return val;
	}
	
	PetscScalar PETScBackend::norm1(const PETScVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_1, &val);
		return val;
	}
	
	PetscScalar PETScBackend::norm_infty(const PETScVector &v) {
		PetscReal val;
		VecNorm(v.implementation(), NORM_INFINITY, &val);
		return val;
	}
	
	PetscScalar PETScBackend::reduce(const PETScVector &vec, const Plus &) {
		//if(!vec.isInitialized()) return 0; // FIXME
		
		PetscScalar result = 0;
		VecSum(vec.implementation(), &result);
		return result;
	}
	
	PetscScalar PETScBackend::reduce(const PETScMatrix &mat, const Plus &op) {
		// assert(false);
		PETScVector rowSum; //FIXME initialize
		Size gs, ls;
		size(mat, gs);
		local_size(mat, ls);
		resize(ls, gs, rowSum);

		MatGetRowSum(mat.implementation(), rowSum.implementation());
		return reduce(rowSum, op);
	}
	
	PetscScalar PETScBackend::reduce(const PETScVector &v, const Min &)
	{
		PetscScalar x;
		VecMin(v.implementation(), nullptr, &x);
		return x;
	}

	PetscScalar PETScBackend::reduce(const PETScMatrix &m, const Min &)
	{
		PetscScalar x;
		MPI_Comm comm = m.communicator();
		
		
		PetscInt grows, gcols, lrows, lcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		MatGetLocalSize(m.implementation(), &lrows, &lcols);
		
		Vec v;
		VecCreate(comm, &v);

		VecSetType(v, VECMPI);
		VecSetSizes(v, lrows, grows);

		MatGetRowMin(m.implementation(), v, nullptr);
		VecMin(v, nullptr, &x);

		VecDestroy(&v);
		return x;
	}

	PetscScalar PETScBackend::reduce(const PETScVector &v, const Max &)
	{
		PetscScalar x;
		VecMax(v.implementation(), nullptr, &x);
		return x;
	}

	PetscScalar PETScBackend::reduce(const PETScMatrix &m, const Max &)
	{
		PetscScalar x;
		PETScVector v;
		VecSetType(v.implementation(), VECMPI);

		PetscInt grows, gcols;
		MatGetSize(m.implementation(), &grows, &gcols);
		VecSetSizes(v.implementation(), PETSC_DECIDE, grows);

		MatGetRowMax(m.implementation(), v.implementation(), nullptr);
		VecMax(v.implementation(), nullptr, &x);
		return x;
	}

	
	// get diagonal of matrix as vector
	bool PETScBackend::diag(PETScVector &vec, const PETScMatrix &mat)
	{
		using std::max;
		PetscInt globalRows, globalColumns;
		
		PetscErrorCode err = MatGetSize(mat.implementation(), &globalRows, &globalColumns);
		PetscInt localRows, localColumns;
		err = MatGetLocalSize(mat.implementation(), &localRows, &localColumns);
		
		PetscInt lenGlobal = max(globalRows, globalColumns);
		PetscInt lenLocal  = max(localRows, localColumns);
		vec.resize(lenLocal, lenGlobal);
		err = MatGetDiagonal(mat.implementation(), vec.implementation());
		return PETScError::Check(err);
	}
	
	bool PETScBackend::diag(PETScSparseMatrix &mat, const PETScVector &vec)
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
		}
		
		
		PetscInt err = MatDiagonalSet( mat.implementation(), vec.implementation(), INSERT_VALUES);
		return PETScError::Check(err);
		
	}

	bool PETScBackend::diag(PETScMatrix &out, const PETScMatrix &in)
	{
		PETScVector vec;
		if(!diag(vec, in)) {return false;}
		return diag(out, vec);
	}

	bool PETScBackend::diag(PETScMatrix &mat, const PETScVector &vec)
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
		return PETScError::Check(err);
		
	}

	bool PETScBackend::mat_diag_shift(PETScMatrix &left, const PetscScalar diag_factor)
	{
		return PETScError::Check(MatShift(left.implementation(), diag_factor));
	}
	
	bool PETScBackend::compare(const Vector &left, const Vector &right, const ApproxEqual &comp) {
		PETScVector diff;
		apply(left, right, Minus(), diff);
		return norm_infty(diff) <= comp.getTol();
	}
	
	bool PETScBackend::compare(const Matrix &left, const Matrix &right, const ApproxEqual &comp) {
		Matrix diff;
		apply(left, right, Minus(), diff);
		return norm2(diff) <= comp.getTol();
	}
	
	bool PETScBackend::apply(const PetscScalar factor, const Vector &vec, const Multiplies &, Vector &result) {
		result = vec;
		VecScale(result.implementation(), factor);
		return true;
	}
	
	bool PETScBackend::apply(const PetscScalar factor, const Matrix &mat, const Multiplies &, Matrix &result) {
		result = mat;
		MatScale(result.implementation(), factor);
		return true;
	}
	
	PetscScalar PETScBackend::dot(const Vector &left, const Vector &right) const {
		PetscScalar result;
		VecDot(left.implementation(), right.implementation(), &result);
		return result;
	}
	
	
	bool PETScBackend::mul(const PETScMatrix &left, const PETScMatrix &right, PETScMatrix &result) {
		
		
		
		Size lsize, rsize;
		size(left, lsize);
		size(right, rsize);
		
		result.initGlobal(lsize.get(0), rsize.get(1));
		
		MatType type;
		MatGetType(right.implementation(), &type);
		
		//            //FIXME: When petsc is fixed for MATMPIDENSE remove branch
		//            if (std::string(MATMPIDENSE) == type) {
		//                Mat sparseleft, sparseright, sparseresult;
		//                MatConvert(left.implementation(), MATMPIAIJ, MAT_INITIAL_MATRIX, &sparseleft);
		//                MatConvert(right.implementation(), MATMPIAIJ, MAT_INITIAL_MATRIX, &sparseright);
		//                if (!PETScError::Check(
		//                        MatMatMult(sparseleft, sparseright, MAT_INITIAL_MATRIX,
		//                                   PETSC_DEFAULT, &sparseresult))) {
		//                    return false;
		//                }
		//
		//                MatConvert(sparseresult, MATMPIDENSE, MAT_INITIAL_MATRIX, &result.implementation());
		//                return true;
		//
		//            }
		
		return PETScError::Check(
								 MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX,
											PETSC_DEFAULT, &result.implementation()));
	}
	
	bool PETScBackend::outer(const Vector &left, const Vector &right, Matrix &result) {
		MPI_Comm comm;
		PetscObjectGetComm((PetscObject)right.implementation(),&comm);
		
		int n_procs = 0, rank = 0;
		MPI_Comm_size(comm, &n_procs);
		MPI_Comm_rank(comm, &rank);
		
		Size lsize, rsize;
		size(left, lsize);
		size(right, rsize);
		const Size result_size = {lsize.get(0), rsize.get(0)};
		
		
		
		const PetscScalar * right_array = nullptr;
		PetscErrorCode err =  VecGetArrayRead(right.implementation(), &right_array);
		
		const PetscScalar * left_array = nullptr;
		err =  VecGetArrayRead(left.implementation(), &left_array);
		
		const Range r_range = range(right);
		const Range l_range = range(left);
		const PetscInt n = rsize.get(0);
		
		std::vector<PetscScalar> recvbuf(n, 0);
		
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
		
		MatAssemblyBegin(result.implementation(), MAT_FINAL_ASSEMBLY);
		for(SizeType i = l_range.begin(); i != l_range.end(); ++i) {
			const PetscScalar l_value = left_array[i-l_range.begin()];
			
			for(SizeType j = 0; j != n; ++j) {
				const PetscScalar r_value = recvbuf.at(j);
				
				MatSetValue(result.implementation(), i, j, l_value*r_value, INSERT_VALUES);
			}
		}
		
		MatAssemblyEnd(result.implementation(), MAT_FINAL_ASSEMBLY);
		VecRestoreArrayRead(left.implementation(), &left_array);
		return true;
	}
	
	void PETScBackend::assignTransposed(PETScMatrix &left, const PETScMatrix &right) {
		if(&left != &right) {
			MatDestroy(&left.implementation());
		} else {
			assert(false);
		}


		MatTranspose(right.implementation(), MAT_INITIAL_MATRIX, &left.implementation());
	}
	
	void PETScBackend::assignToRange(PETScMatrix & /*left*/, const PETScMatrix &/*right*/, const Range &/*globalRowRange*/,
									 const Range &/*globalColRange*/)
	{
		
		assert(false); //TODO
		
	}
	
	void PETScBackend::assignToRange( PETScMatrix &/*left*/, const Identity &/**/, const Range &/*globalRowRange*/,
									 const Range &/*globalColRange*/) {
		assert(false); //TODO
	}
	
	void PETScBackend::assignToRange( PETScVector &/*left*/, const PETScVector &/*right*/, const Range &/*globalRowRange*/,
									 const Range &/*globalColRange*/) {
		assert(false); //TODO
	}
	
	bool PETScBackend::vec2mat(const Vector & v, Matrix & m , const bool  transpose ) {
		
		PetscInt n;
		VecGetSize(v.implementation(), &n);
		Range r = range(v);
		
		Matrix mtr;
		
		Matrix &work = transpose ? mtr : m;
		this->build(work, Size({n, 1}), Values<PetscScalar>(0));
		
		MatAssemblyBegin(work.implementation(), MAT_FINAL_ASSEMBLY);
		
		PetscScalar zero = 0;
		for (PetscInt i = r.begin(); i < r.end(); ++i) {
			PetscScalar val = get(v, i);
			set(work, i, zero, val);
		}
		
		MatAssemblyEnd(work.implementation(), MAT_FINAL_ASSEMBLY);
		
		if (transpose) {
			MatTranspose(mtr.implementation(), MAT_INITIAL_MATRIX, &m.implementation());
		}
		
		return true;
	}
	
	bool PETScBackend::gemv(  const PetscScalar /*alpha */, const PETScMatrix &/*left*/, const PETScVector &/*right*/,
							const PetscScalar/* beta*/, PETScVector &/*result*/) {
		assert(false); //TODO
		return false;
	}
	
	
	// this is not very correct yet ...
	bool PETScBackend::build_local_diag_block(PETScSerialSparseMatrix &left, const PETScSparseMatrix &right)
	{
		
		Mat M;
		MatGetDiagonalBlock(right.implementation(), &M);
		convert(M, left);
		
		return true;
	}
	
	bool PETScBackend::triple_product_PtAP(const PETScMatrix & A, const PETScMatrix & P, PETScMatrix & result)
	{
		if(&result.implementation() != &A.implementation() && &result.implementation() != &P.implementation()) {
			MatDestroy(&result.implementation());
		} else {
			std::cerr << "[Warning] not handled case in triple_product_PtAP" << std::endl;
		}

		MatPtAP(A.implementation(), P.implementation(), MAT_INITIAL_MATRIX, 1.0, &result.implementation()); 
		return true; 
	}
	


	bool PETScBackend::triple_product(const PETScMatrix & A, const PETScMatrix & B, const PETScMatrix & C, PETScMatrix & result)
	{
		if(&result.implementation() != &A.implementation() && &result.implementation() != &B.implementation() && &result.implementation() != &C.implementation()) {
			MatDestroy(&result.implementation());
		}
		
		MatMatMatMult(A.implementation(), B.implementation(), C.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()); 
		return true; 
	}



	bool PETScBackend::is_nan_or_inf(const PETScVector & X)
	{
		PetscInt m; 
        const PetscScalar *x;
        VecGetLocalSize(X.implementation(), &m);
        VecGetArrayRead(X.implementation(), &x);

        PetscErrorCode error; 

        for (PetscInt i=0; i<m; i++) 
        {
            error  =PetscIsInfOrNanScalar(x[i]); 
            if(error ==1)
            	return 1;
        }
        VecRestoreArrayRead(X.implementation(), &x);
		return 0; 
	}

	// redistribution of local sizes of vector
	// TODO: can be done also for matrices 
	//       can be done also based on local sizes, no provided vector
	bool PETScBackend::build_local_redistribute(const PETScVector &x_from, const PETScVector &shape_vec, PETScVector &result)
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
        
        return true; 		
	}


	bool PETScBackend::aux_zaxpy(const PetscScalar scaleFactor, const PETScVector &left,
								 PETScVector &result) {
		
		VecAXPY(result.implementation(), scaleFactor, left.implementation());
		return true;
	}

	bool PETScBackend::waxpby(const PetscScalar a, const Vector &x, const PetscScalar &b, const Vector &y, Vector &result)
	{
		if (result.implementation() == y.implementation()) {
			VecAXPBY(result.implementation(), a, b, x.implementation());
		} else {
			if(result.implementation() == x.implementation()) {
				VecAXPBY(result.implementation(), b, a, y.implementation());
			} else {
				result = y;
				VecAXPBY(result.implementation(), a, b, x.implementation());
			}
		}
		
		return true;
	}

	
	bool PETScBackend::aux_zaxpy(const PetscScalar scaleFactor, const PETScMatrix &left,
								 PETScMatrix &result) {
		MatAXPY(result.implementation(), scaleFactor, left.implementation(), DIFFERENT_NONZERO_PATTERN);
		return true;
	}

	bool PETScBackend::diag_scale_right(const Matrix &m, const Vector &diag, Matrix &result)
	{
		assign(result, m);
		return PETScError::Check( MatDiagonalScale(result.implementation(), nullptr, diag.implementation()) );
	}

	bool PETScBackend::diag_scale_left(const Vector &diag, const Matrix &m, Matrix &result)
	{
		assign(result, m);
		return PETScError::Check( MatDiagonalScale(result.implementation(), diag.implementation(), nullptr) );
	}

	bool PETScBackend::diag_scale_left(const Vector &diag, const Vector &m,  Vector &result)
	{
        Size gs, ls;
        size(m, gs);
        local_size(m, ls);
        resize(ls, gs, result);
        return PETScError::Check( VecPointwiseMult(diag.implementation(), m.implementation(), result.implementation()) );
	}


	PetscScalar PETScBackend::trace(const Matrix &mat)
	{
		PetscScalar ret;
		MatGetTrace(mat.implementation(), &ret);
		return ret;
	}

	bool PETScBackend::apply_tensor_reduce(const Matrix &mat, const Plus &, const int dim, Vector &result)
	{
		PETScVector rowSum; //FIXME initialize
		Size gs, ls;
		size(mat, gs);
		local_size(mat, ls);
		resize(ls, gs, result);

		if(dim == 1) {
			MatGetRowSum(mat.implementation(), result.implementation());
		} else {
			//FIXME implement own
			assert(false && "not available in pestsc");
			return false;
		}

		return true;
	}

	bool PETScBackend::apply_tensor_reduce(const Matrix &mat, const Min &, const int dim, Vector &result)
	{
		// PetscScalar x;
		PetscInt grows, gcols;
		MatGetSize(mat.implementation(), &grows, &gcols);

		if (dim == 1) {
			VecDestroy(&result.implementation());
			VecCreateMPI(mat.communicator(), PETSC_DECIDE, grows, &result.implementation());
			MatGetRowMin(mat.implementation(), result.implementation(), nullptr);
		} else {
			//FIXME implement own
			// VecDestroy(&result.implementation());
			// VecCreateMPI(mat.communicator(), PETSC_DECIDE, gcols, &result.implementation());
			assert(false && "not available in petsc");
			return false;
		}

		return true;
	}

	bool PETScBackend::apply_tensor_reduce(const Matrix &mat, const Max &, const int dim, Vector &result)
	{
		// PetscScalar x;
		PetscInt grows, gcols;
		MatGetSize(mat.implementation(), &grows, &gcols);

		if (dim == 1) {
			VecDestroy(&result.implementation());
			VecCreateMPI(mat.communicator(), PETSC_DECIDE, grows, &result.implementation());
			MatGetRowMax(mat.implementation(), result.implementation(), nullptr);
		} else {
			//FIXME implement own
			// VecDestroy(&result.implementation());
			// VecCreateMPI(mat.communicator(), PETSC_DECIDE, gcols, &result.implementation());
			assert(false && "not available in petsc");
			return false;
		}

		return true;
	}

	bool PETScBackend::inverse(const Matrix &mat, Matrix &result)
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
		if(!PETScError::Check( MatMatSolve(L.implementation(), I.implementation(), result.implementation()) )) {
			return false;
		}

		// MatAssemblyBegin(result.implementation(), MAT_FINAL_ASSEMBLY);
		// MatAssemblyEnd(result.implementation(), MAT_FINAL_ASSEMBLY);


		ISDestroy(&isr);
		ISDestroy(&isc);
		return true;
	}


}

