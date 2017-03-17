#ifndef UTOPIA_UTOPIA_PETSCBACKEND_HPP
#define UTOPIA_UTOPIA_PETSCBACKEND_HPP

#include "utopia_PETScTraits.hpp"
#include "utopia_Core.hpp"
#include "utopia_Factory.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_Base.hpp"
#include "utopia_ScalarBackend.hpp"



namespace utopia
{
	class CompatibleMatPair {
	public:
		CompatibleMatPair(const MPI_Comm comm, const Mat &left, const Mat &right);
		inline const Mat &left() const
		{
			return left_;
		}

		inline const Mat &right() const
		{
			return right_;
		}

		~CompatibleMatPair();

	private:
		static bool is_sparse(const Mat &mat);

		Mat left_;
		Mat right_;
		bool must_destroy_left_;
		bool must_destroy_right_;
	};

	bool is_matlab_file(const std::string &path);


	class PETScBackend : public ScalarBackend<PetscScalar>
	{

	public:
		typedef PetscScalar Scalar;
		typedef PETScVector Vector;
		typedef PETScMatrix Matrix;


		using ScalarBackend<PetscScalar>::apply;
		using ScalarBackend<PetscScalar>::zaxpy;

		template<class LorRValueMatrix>
		void assign(PETScMatrix &left, LorRValueMatrix &&right)
		{
			left = std::forward<LorRValueMatrix>(right);
		}

		template<class LorRValueMatrix>
		void assign(PETScSparseMatrix &left, LorRValueMatrix &&right)
		{
			left = std::forward<LorRValueMatrix>(right);
		}

		template<class LorRValueVector>
		void assign(PETScVector &left, LorRValueVector &&right)
		{
			left = std::forward<LorRValueVector>(right);
		}

		void clear(PETScMatrix &mat) const;

		//void structure(Vec vec, Structure &structure)

		bool convert(Vec vec, PETScVector &wrapper);
		bool convert(Mat mat, PETScMatrix &wrapper);
		bool convert(Mat mat, PETScSparseMatrix &wrapper);
		bool convert(PETScVector &wrapper, Vec vec);

		template<int FillType>
		bool convert(PETScGenericMatrix<FillType> &wrapper, Mat mat)
		{
			MatCopy(wrapper.implementation(), mat, SAME_NONZERO_PATTERN);
			return true;
		}

		bool wrap(Mat mat, PETScSparseMatrix &wrapper);
		bool wrap(Vec vec, PETScVector &wrapper);

		Range range(const PETScVector &v);

		template<int FillType>
		Range rowRange(const PETScGenericMatrix<FillType> &m) {
			PetscInt rbegin, rend;
			MatGetOwnershipRange(m.implementation(), &rbegin, &rend);
			assert(Range(rbegin, rend).valid());
			return Range(rbegin, rend);
		}

		template<int FillType>
		Range colRange(const PETScGenericMatrix<FillType> &m) {
			PetscInt grows, gcols;
			MatGetSize(m.implementation(), &grows, &gcols);
			return Range(0, gcols);
		}

		template<int FillType>
		bool size(const PETScGenericMatrix<FillType> &m, Size &size) const {
			PetscInt grows, gcols;
			size.setDims(2);
			MatGetSize(m.implementation(), &grows, &gcols);
			size.set(0, grows);
			size.set(1, gcols);
			return true;
		}

		template<int FillType>
		bool local_size(const PETScGenericMatrix<FillType> &mat, Size &size) const {
			PetscInt n, m;
			size.setDims(2);
			MatGetLocalSize(mat.implementation(), &n, &m);
			size.set(0, n);
			size.set(1, m);
			return true;
		}

		bool size(const PETScVector &m, Size &size) const;
		bool local_size(const PETScVector &m, Size &size) const;

		template<int FillType>
		void resize(const Size &s, PETScGenericMatrix<FillType> &mat)
		{
			PetscInt r, c;
			MatGetSize(mat.implementation(), &r, &c);
			if(r == s.get(0) && c == s.get(1)) {
				return;
			}

			MatType type;
			PETScError::Check( MatGetType(mat.implementation(), &type) );

			mat.init();
			MatSetSizes(mat.implementation(), PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1));
			MatSetType(mat.implementation(), type);
		}

		template<int FillType>
		void resize(const Size &local_s, const Size &global_s, PETScGenericMatrix<FillType> &mat)
		{
			PetscInt r, c, R, C;
			MatGetLocalSize(mat.implementation(), &r, &c);
			MatGetSize(mat.implementation(), &R, &C);
			if(r == local_s.get(0) && c == local_s.get(1) && R == global_s.get(0) && C == global_s.get(1)) {
				return;
			}

			MatType type;
			PETScError::Check( MatGetType(mat.implementation(), &type) );
			mat.init();
			MatSetSizes(mat.implementation(), local_s.get(0), local_s.get(1), global_s.get(0), global_s.get(1));
			MatSetType(mat.implementation(), type);
		}

		void resize(const Size &size, PETScVector &mat);
		void resize(const Size &s_local, const Size &s_global, PETScVector &vec);

		template<int FillTypeLeft, int FillTypeRight>
		void assignFromRange(PETScGenericMatrix<FillTypeLeft> &left, const PETScGenericMatrix<FillTypeRight> &right,
				const Range &globalRowRange, const Range &globalColRange)
		{
			assert(!globalRowRange.empty());

			MPI_Comm comm = right.communicator();
			left.setCommunicator(comm);
			left.init({PETSC_DECIDE, PETSC_DECIDE}, {globalRowRange.extent(), globalColRange.extent()});


			Range rr = rowRange(right).intersect(globalRowRange);
			Range cr = colRange(right).intersect(globalColRange);

			// rr might be invalid !!! so use < to compare with range
			//            if(!rr.valid()) std::cout << rr << std::endl;


			//use MAT_FLUSH_ASSEMBLY to not make petsc crash if nothing is added
			MatAssemblyBegin(left.implementation(), MAT_FLUSH_ASSEMBLY);

			PetscInt r = 0, c = 0;
			for (PetscInt rIt = rr.begin(); rIt < rr.end(); ++rIt) {
			 r = rIt - globalRowRange.begin();
			 for (PetscInt cIt = cr.begin(); cIt < cr.end(); ++cIt) {
				 c = cIt - globalColRange.begin();
				 set(left, r, c, get(right, rIt, cIt));
			 }
			}

			MatAssemblyEnd(left.implementation(), MAT_FLUSH_ASSEMBLY);

			//To finalize the matrix
			MatAssemblyBegin(left.implementation(), MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(left.implementation(), MAT_FINAL_ASSEMBLY);

		}

		void assignFromRange(PETScVector &left, const PETScVector &right, const Range &globalRowRange,
							 const Range & /*globalColRange */);
		// read matrix
		template<int FillType>
		bool read(const std::string &path, PETScGenericMatrix<FillType> &Mat_A)
		{
			Mat_A.init();
			PetscViewer fd;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_READ, &fd);
			bool status;
			status =  PETScError::Check( MatLoad(Mat_A.implementation(),fd) );
			PetscViewerDestroy(&fd);
			return status;
		}

		// write matrix
		template<int FillType>
		bool write(const std::string &path, const PETScGenericMatrix<FillType> &Mat_A)
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
		bool write(const std::string &path, const PETScVector &Vec_A);

		// monitor for cyrill
		bool monitor(const long & it, PETScVector &Vec_A);

		template<int FillType>
		bool monitor(const long & it, PETScGenericMatrix<FillType> &Mat_A)
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

		template<int FillType>
		PetscScalar get_global_nnz(PETScGenericMatrix<FillType> &Mat_A) {
			MatInfo info;
			MatGetInfo(Mat_A.implementation(), MAT_GLOBAL_SUM, &info);
			return info.nz_used;
		}

		template<int FillType>
		PetscScalar get_local_nnz(PETScGenericMatrix<FillType> &Mat_A) {
			MatInfo        info;
			MatGetInfo(Mat_A.implementation(), MAT_LOCAL, &info);
			return info.nz_used;
		}

		// call to MatZeroRows
		template<int FillType>
		bool set_zero_rows(PETScGenericMatrix<FillType> &Mat_A, const std::vector<int> &index) {
			return PETScError::Check(MatZeroRows(Mat_A.implementation(), index.size(), &index[0], 0, NULL, NULL));
		}

		template<int FillType>
		bool apply_BC_to_system(PETScGenericMatrix<FillType> & A, PETScVector& x, PETScVector& rhs, const std::vector<int> &index) {
			return PETScError::Check(MatZeroRows(A.implementation(), index.size(), &index[0], 1.0, x.implementation(), rhs.implementation()));
		}

		// read vector
		bool read(const std::string &path, PETScVector &Vec_A);

		template<int FillType>
		bool scal(const PetscScalar scaleFactor, const PETScGenericMatrix<FillType> &left, PETScGenericMatrix<FillType> &right) {
			if(&left == &right) {
				MatScale(right.implementation(), scaleFactor);
			} else {
				assert(false); // FIXME
				return false;
			}

			return true;
		}
		// bool apply(const PetscScalar scaleFactor, const PETScMatrix &m, const Multiplies &, PETScMatrix &result)
		// {
		// 	return scal(scaleFactor, m, result);
		// }

		template<class Tensor, class LorRValueTensor, class ResultTensor>
		bool zaxpy(const PetscScalar scaleFactor, const Tensor &left, LorRValueTensor &&right, ResultTensor &result)
		{
			result = std::forward<LorRValueTensor>(right);

			if (result.implementation() == left.implementation()) {
				//TODO
				std::cerr << "axpy" << std::endl;
			}

			return aux_zaxpy(scaleFactor, left, result);
		}

		bool aux_zaxpy(const PetscScalar scaleFactor, const PETScVector &left, PETScVector &result);

		template<int FillType>
		bool aux_zaxpy(const PetscScalar scaleFactor, const PETScGenericMatrix<FillType> &left,
									 PETScGenericMatrix<FillType> &result) {
			MatAXPY(result.implementation(), scaleFactor, left.implementation(), DIFFERENT_NONZERO_PATTERN);
			return true;
		}

		template<int FillType>
		bool gemm(const PetscScalar alpha, const PETScGenericMatrix<FillType> &left, const PETScGenericMatrix<FillType> &right,
				bool transpose_left, bool transpose_right, const PetscScalar beta, PETScGenericMatrix<FillType> &result) {
			//FIXME only works for beta == 0 for the moment
			assert(fabs(beta) < 1e-16);
			//FIXME only works for alpha == 1 for the moment
			assert(fabs(alpha - 1) < 1e-16);

			CompatibleMatPair mat_pair(left.communicator(), left.implementation(), right.implementation());
			auto l = mat_pair.left();
			auto r = mat_pair.right();

			result.init();
			bool ok = false;
			if(transpose_left && !transpose_right) {
				ok = PETScError::Check(MatTransposeMatMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()));
			} else if(!transpose_left && transpose_right) {
				ok = PETScError::Check(MatMatTransposeMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()));
			} else if(!transpose_left && !transpose_right) {
				ok = PETScError::Check(MatMatMult(l, r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()));
			} else {
				assert(transpose_left && transpose_right);
				PETScGenericMatrix<FillType> temp;

				if(!PETScError::Check( MatMatMult(r, l, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation()) )) {
					ok = false;
				} else {
					ok = transpose(temp, result);
				}
			}

			return ok;
		}

		bool gemm(const PetscScalar /*alpha*/, const PETScMatrix &/*left*/, const PETScMatrix &/*right*/, const PetscScalar /* beta*/,
				  PETScMatrix & /*result*/);

		bool gemm(const PetscScalar /*alpha */, const PETScVector &/*left*/, const PETScMatrix &/*right*/,
				  bool /*transpose_left*/, bool /*transpose_right*/, const PetscScalar /*beta*/, PETScMatrix &/*result*/);

		template<int FillType>
		bool gemm(const double scaleFactor, const PETScGenericMatrix<FillType> &left, const Vector &right,
				bool transpose_left,  bool transpose_right, const double beta, Vector &result)
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


		template<class Left, class Right, class Result>
		bool gem(const double /*scaleFactor*/, const Left &/*left*/, const Right &/*right*/, bool /*transpose_left*/,
				 bool /*transpose_right*/, const double /*beta*/, Result &/*result*/);


		template<int FillType>
		bool transpose(const PETScGenericMatrix<FillType> &mat, PETScGenericMatrix<FillType> &result)
		{
			if(&mat != &result) {
				result.init();
			} else {
				assert(false);
			}

			return PETScError::Check(MatTranspose(mat.implementation(), MAT_INITIAL_MATRIX, &result.implementation()));
		}


		bool scal(const PetscScalar /*scale_factor*/, const PETScVector &/*v*/, PETScVector &/*result */);

		template<class Tensor>
		void build(Tensor &t, const Size &s, const Resize &)
		{
			//FIXME switch with more efficient version without zeros
			build(t, s, Zeros());
			// resize(s, t);
		}

		// Too different to make a template for these
		void build(PETScMatrix &m, const Size &size, const Identity &);
		void build(PETScSparseMatrix &m, const Size &size, const Identity &);
		void build(PETScMatrix &m, const Size &size, const LocalIdentity &);
		void build(PETScSparseMatrix &m, const Size &size, const LocalIdentity &);

		// TODO - find a way to pass NNZ to PETScSparseMatrix (or its Allocator)
		void build(PETScSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz);

		void build(PETScSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz);
		void build(PETScSparseMatrix &m, const Size &size, const LocalRowNNZ<PetscInt> &nnz);
		/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
		void build(PETScMatrix  &m, const Size &size, const LocalNNZ<PetscInt> & /*nnz */);

		/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
		void build(PETScMatrix  &m, const Size &size, const NNZ<PetscInt> &/*nnz*/);

		template<int FillType>
		void build(PETScGenericMatrix<FillType> &m, const Size &size, const Zeros &) {
			build(m, size, Values<PetscScalar>(0));
		}

		template<int FillType>
		void build(PETScGenericMatrix<FillType> &m, const Size &size, const LocalZeros &) {
			build(m, size, LocalValues<PetscScalar>(0));
		}

		void build(PETScVector &v, const Size &size, const Zeros &);
		void build(PETScVector &v, const Size &size, const LocalZeros &);

		template<int FillType>
		void build(PETScGenericMatrix<FillType> &m, const Size &size, const Values<PetscScalar> &values) {
			m.init({PETSC_DECIDE, PETSC_DECIDE}, size);

			PetscInt rbegin, rend;
			MatGetOwnershipRange(m.implementation(), &rbegin, &rend);

			PetscInt grows, gcols;
			MatGetSize(m.implementation(), &grows, &gcols);

			if (FillType == FillType::SPARSE)
				if (MatSeqAIJSetPreallocation(m.implementation(), gcols, NULL)) assert(false);

			MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);

			const PetscScalar v = values.value();
			for (PetscInt i = rbegin; i < rend; ++i) {
				for (PetscInt j = 0; j < gcols; ++j) {
					MatSetValue(m.implementation(), i, j, v, INSERT_VALUES);
				}
			}

			MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
		}

		void build(PETScVector &v, const Size &local_size, const Size &&global_size, const Values<PetscScalar> &values);
		void build(PETScVector &v, const Size &size, const Values<PetscScalar> &values);

		template<int FillType>
		void build(PETScGenericMatrix<FillType> &m, const Size &size, const LocalValues<PetscScalar> &values) {
			m.init(size, {PETSC_DETERMINE, PETSC_DETERMINE});

			//TODO check if it can be simplified using known information.

			PetscInt rbegin, rend;
			MatGetOwnershipRange(m.implementation(), &rbegin, &rend);

			PetscInt grows, gcols;
			MatGetSize(m.implementation(), &grows, &gcols);

			MatAssemblyBegin(m.implementation(), MAT_FINAL_ASSEMBLY);

			if (FillType == FillType::SPARSE)
				if (MatSeqAIJSetPreallocation(m.implementation(), gcols, NULL)) assert(false);

			const PetscScalar v = values.value();
			for (PetscInt i = rbegin; i < rend; ++i) {
				for (PetscInt j = 0; j < gcols; ++j) {
					MatSetValue(m.implementation(), i, j, v, ADD_VALUES);
				}
			}

			MatAssemblyEnd(m.implementation(), MAT_FINAL_ASSEMBLY);
		}

		void build(PETScVector &v, const Size &size, const LocalValues<PetscScalar> &values);

		void set(PETScVector &v, const PetscInt index, PetscScalar value);
		void add(PETScVector &v, const PetscInt index, PetscScalar value);

		void set(PETScVector &v, const std::vector<PetscInt> indices, const std::vector<PetscScalar> values);

		template<int FillType>
		void set(PETScGenericMatrix<FillType> &v, const PetscInt row, const PetscInt col, PetscScalar value) {
			MatSetValues(v.implementation(), 1, &row, 1, &col, &value, INSERT_VALUES);
		}

		template<int FillType>
		void add(PETScGenericMatrix<FillType> &m, const PetscInt row, const PetscInt col, PetscScalar value) {
			MatSetValues(m.implementation(), 1, &row, 1, &col, &value, ADD_VALUES);
		}


		template<class Tensor>
		void readLock(const Tensor &) {}

		template<class Tensor>
		void readUnlock(const Tensor &) {}

		void writeLock(PETScVector &vec);
		void writeUnlock(PETScVector &vec);

		template<int FillType>
		void writeLock(const PETScGenericMatrix<FillType> &mat) {
			MatAssemblyBegin(mat.implementation(), MAT_FINAL_ASSEMBLY);
		}

		template<int FillType>
		void writeUnlock(const PETScGenericMatrix<FillType> &mat) {
			MatAssemblyEnd(mat.implementation(), MAT_FINAL_ASSEMBLY);
		}

		template<int FillType>
		void set(PETScGenericMatrix<FillType> &v, const std::vector<PetscInt> rows, const std::vector<PetscInt> cols,
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

		PetscScalar get(const PETScVector &v, const PetscInt index);

		template<int FillType>
		PetscScalar get(const PETScGenericMatrix<FillType> &v, const PetscInt row, const PetscInt col) {
			PetscScalar value;
			MatGetValues(v.implementation(), 1, &row, 1, &col, &value);
			return value;
		}

		template<int FillType>
		bool apply(const PETScGenericMatrix<FillType> &left, const PETScVector &right, const Multiplies &, PETScVector &result) {
			PetscInt grows, gcols;
			MatGetSize(left.implementation(), &grows, &gcols);

			PetscInt rows, cols;
			MatGetLocalSize(left.implementation(), &rows, &cols);

			result.setCommunicator(left.communicator());

			result.init({rows}, {grows});
			VecAssemblyBegin(result.implementation());
			VecAssemblyEnd(result.implementation());

			MatMult(left.implementation(), right.implementation(), result.implementation());

			return true;
		}

		template<int FillTypeLeft, int FillTypeRight, int FillTypeResult>
		bool apply(const PETScGenericMatrix<FillTypeLeft> &left, const PETScGenericMatrix<FillTypeRight> &right, const Multiplies &,
				PETScGenericMatrix<FillTypeResult> &result)
		{
			if(&right.implementation() != &result.implementation() || &left.implementation() !=  &result.implementation()) {
				result.init();
			} else {
				assert(false);
			}

			MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
			return true;
		}

		template<int FillType>
		bool mat_mult_add(const PETScGenericMatrix<FillType> &m, const PETScVector &right, const PETScVector &left, PETScVector &result) {
			if (&right.implementation() == &result.implementation() || &left.implementation() ==  &result.implementation()) {
				assert(false);
			}

			PetscInt size, gsize;
			VecGetSize(right.implementation(), &gsize);
			VecGetLocalSize(right.implementation(), &size);

			result.setCommunicator(right.communicator());

			result.init({size}, {gsize});
			VecAssemblyBegin(result.implementation());
			VecAssemblyEnd(result.implementation());

			MatMultAdd(m.implementation(), right.implementation(), left.implementation(), result.implementation());
			return true;
		}

		template<int FillType>
		bool mat_multT_add(const PETScGenericMatrix<FillType> &m, const PETScVector &right, const PETScVector &left, PETScVector &result) {
			if (&right.implementation() == &result.implementation() || &left.implementation() ==  &result.implementation()) {
				assert(false);
			}

			PetscInt size, gsize;
			VecGetSize(right.implementation(), &gsize);
			VecGetLocalSize(right.implementation(), &size);

			result.setCommunicator(right.communicator());

			result.init({size}, {gsize});
			VecAssemblyBegin(result.implementation());
			VecAssemblyEnd(result.implementation());

			MatMultTransposeAdd(m.implementation(), right.implementation(), left.implementation(), result.implementation());
			return true;
		}

		bool apply(const PETScVector &left, const PETScVector &right, const EMultiplies &, PETScVector &result);
		bool apply(const PETScVector &left, const PETScVector &right, const Divides &, PETScVector &result);

		// reciprocal
		bool apply(const PETScVector &vec, const Reciprocal<double> &reciprocal, PETScVector &result);
		PetscScalar norm2(const PETScVector &v);

		template<int FillType>
		PetscScalar norm2(const PETScGenericMatrix<FillType> &m) {
			PetscReal val;
			MatNorm(m.implementation(), NORM_FROBENIUS, &val);
			return val;
		}

		PetscScalar norm1(const PETScVector &v);
		PetscScalar norm_infty(const PETScVector &v);
		PetscScalar reduce(const PETScVector &vec, const Plus &);

		template<int FillType>
		PetscScalar reduce(const PETScGenericMatrix<FillType> &mat, const Plus &op) {
			PETScVector rowSum;
			Size gs, ls;
			size(mat, gs);
			local_size(mat, ls);
			resize(ls, gs, rowSum);

			MatGetRowSum(mat.implementation(), rowSum.implementation());
			return reduce(rowSum, op);
		}

		PetscScalar reduce(const PETScVector &, const Min &);
		PetscScalar reduce(const PETScVector &, const Max &);

		template<int FillType>
		PetscScalar reduce(const PETScGenericMatrix<FillType> &m, const Min &) {
			PetscScalar x;
			PETScVector v;
			VecSetType(v.implementation(), VECMPI);

			PetscInt grows, gcols;
			MatGetSize(m.implementation(), &grows, &gcols);
			VecSetSizes(v.implementation(), PETSC_DECIDE, grows);

			MatGetRowMin(m.implementation(), v.implementation(), nullptr);
			VecMin(v.implementation(), nullptr, &x);
			return x;
		}

		template<int FillType>
		PetscScalar reduce(const PETScGenericMatrix<FillType> &m, const Max &) {
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
		template<int FillType>
		bool diag(PETScVector &vec, const PETScGenericMatrix<FillType> &mat) {
			using std::max;
			PetscInt globalRows, globalColumns;

			PetscErrorCode err = MatGetSize(mat.implementation(), &globalRows, &globalColumns);
			PetscInt localRows, localColumns;
			err = MatGetLocalSize(mat.implementation(), &localRows, &localColumns);

			PetscInt lenGlobal = max(globalRows, globalColumns);
			PetscInt lenLocal  = max(localRows, localColumns);
			vec.init({lenLocal}, {lenGlobal});
			err = MatGetDiagonal(mat.implementation(), vec.implementation());
			return PETScError::Check(err);
		}

		template<int FillType>
		bool diag(PETScGenericMatrix<FillType> &mat, const PETScVector &vec)
		{
			// I do not think, this is needed
			// because, doesnt run in parallel properly
			// u can not change size of matrix after it was initialized...
			// or at least it seems like ...
			//!!! FIXME needs to check if matrix has been allocated already

			mat.init();

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

		template<int FillTypeLeft, int FillTypeRight>
		bool diag(PETScGenericMatrix<FillTypeLeft> &out, const PETScGenericMatrix<FillTypeRight> &in) {
			PETScVector vec;
			if(!diag(vec, in)) {return false;}
			return diag(out, vec);
		}

		template<int FillType>
		bool mat_diag_shift(PETScGenericMatrix<FillType> &left, const PetscScalar diag_factor) {
			return PETScError::Check(MatShift(left.implementation(), diag_factor));
		}

		bool compare(const Vector &left, const Vector &right, const ApproxEqual &comp);

		template<int FillType>
		bool compare(const PETScGenericMatrix<FillType> &left, const PETScGenericMatrix<FillType> &right, const ApproxEqual &comp) {
			PETScGenericMatrix<FillType> diff;
			apply(left, right, Minus(), diff);
			return norm2(diff) <= comp.getTol();
		}


		bool apply(const PetscScalar factor, const Vector &vec, const Multiplies &, Vector &result);

		template<int FillType>
		bool apply(const PetscScalar factor, const PETScGenericMatrix<FillType> &mat, const Multiplies &, PETScGenericMatrix<FillType> &result) {
			result = mat;
			MatScale(result.implementation(), factor);
			return true;
		}

		bool apply(const Vector &vec, const Abs &, Vector &result);
		// bool apply(const Matrix &mat, const Abs &, Matrix &result);

		PetscScalar dot(const Vector &left, const Vector &right) const;


		template<int FillTypeLeft, int FillTypeRight, int FillTypeResult>
		bool mul(const PETScGenericMatrix<FillTypeLeft> &left, const PETScGenericMatrix<FillTypeRight> &right, PETScGenericMatrix<FillTypeResult> &result) {
			Size lsize, rsize;
			size(left, lsize);
			size(right, rsize);

			MatSetSizes(result.implementation(), PETSC_DECIDE, PETSC_DECIDE, lsize.get(0), rsize.get(1));

			MatType type;
			MatGetType(right.implementation(), &type);

			return PETScError::Check(
				MatMatMult(left.implementation(), right.implementation(), MAT_INITIAL_MATRIX,
					PETSC_DEFAULT, &result.implementation()));
		}

		bool outer(const Vector &left, const Vector &right, Matrix &result);

		template<int FillType>
		void assignTransposed(PETScGenericMatrix<FillType> &left, const PETScGenericMatrix<FillType> &right) {
			if(&left != &right) {
				left.init();
			} else {
				assert(false);
			}

			MatTranspose(right.implementation(), MAT_INITIAL_MATRIX, &left.implementation());
		}

		void assignToRange(PETScMatrix & /*left*/, const PETScMatrix &/*right*/, const Range &/*globalRowRange*/,
						   const Range &/*globalColRange*/);

		void assignToRange( PETScMatrix &/*left*/, const Identity &/**/, const Range &/*globalRowRange*/,
						   const Range &/*globalColRange*/);

		void assignToRange( PETScVector &/*left*/, const PETScVector &/*right*/, const Range &/*globalRowRange*/,
						   const Range &/*globalColRange*/);

		bool vec2mat(const Vector & v, Matrix & m , const bool  transpose );
		bool gemv(  const PetscScalar /*alpha */, const PETScMatrix &/*left*/, const PETScVector &/*right*/,
				  const PetscScalar/* beta*/, PETScVector &/*result*/);

		// this is not very correct yet ...
		bool build_local_diag_block(PETScSerialSparseMatrix &left, const PETScSparseMatrix &right);
		bool build_local_redistribute(const PETScVector &, const PETScVector &, PETScVector &);

		template<int FillType>
		bool triple_product_PtAP(const PETScGenericMatrix<FillType> &A, const PETScGenericMatrix<FillType> &P, PETScGenericMatrix<FillType> &result) {
			if(&result.implementation() != &A.implementation() && &result.implementation() != &P.implementation()) {
				result.init();
			} //else FIXME

			MatPtAP(A.implementation(), P.implementation(), MAT_INITIAL_MATRIX, 1.0, &result.implementation());
			return true;
		}

		template<int FillType>
		bool triple_product(const PETScGenericMatrix<FillType> & A, const PETScGenericMatrix<FillType> & B, const PETScGenericMatrix<FillType> & C, PETScGenericMatrix<FillType> & result) {
			if(&result.implementation() != &A.implementation() && &result.implementation() != &B.implementation() && &result.implementation() != &C.implementation()) {
				result.init();
			}

			MatMatMatMult(A.implementation(), B.implementation(), C.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
			return true;
		}

		bool is_nan_or_inf(const PETScVector &);


		template<class Tensor>
		void readAndWriteLock(const Tensor &) {}

		template<class Tensor>
		void readAndWriteUnlock(const Tensor &) {}

		template<class LeftTensor, class RightTensor, class ResultTensor>
		inline bool apply(LeftTensor &&left, RightTensor &&right, const Plus &, ResultTensor &result) {
			return zaxpy(1.0, std::forward<LeftTensor>(left), std::forward<RightTensor>(right), result);
		}

		template<class LeftTensor, class RightTensor, class ResultTensor>
		inline bool apply(LeftTensor &&left, RightTensor &&right, const Minus &, ResultTensor &result) {
			return zaxpy(-1.0, std::forward<RightTensor>(right), std::forward<LeftTensor>(left), result);
		}

		template<int FillType>
		bool diag_scale_right(const PETScGenericMatrix<FillType> &m, const Vector &diag, PETScGenericMatrix<FillType> &result) {
			assign(result, m);
			return PETScError::Check( MatDiagonalScale(result.implementation(), nullptr, diag.implementation()) );
		}

		template<int FillType>
		bool diag_scale_left(const Vector &diag, const PETScGenericMatrix<FillType> &m, PETScGenericMatrix<FillType> &result) {
			assign(result, m);
			return PETScError::Check( MatDiagonalScale(result.implementation(), diag.implementation(), nullptr) );
		}

		bool diag_scale_left(const Vector &diag, const Vector &m, Vector &result);

		template<class Operation>
		bool apply(const Vector &v, const Operation &, Vector &result)
		{
			Range r = range(v);

			Size gs, ls;
			size(v, gs);
			local_size(v, ls);
			VecSetType(result.implementation(), VECMPI);
			VecSetSizes(result.implementation(), ls.get(0), gs.get(0));

			writeLock(result);
			readLock(v);

			auto fun = Operation:: template Fun<PetscScalar>();

			for(PetscInt i = r.begin(); i < r.end(); ++i) {
				PetscScalar value = get(v, i);
				VecSetValue(result.implementation(), i, fun(value), INSERT_VALUES);
			}

			writeUnlock(result);
			readUnlock(v);
			return true;
		}


		bool apply(const Vector &v, const Minus &, Vector &result)
		{
			Range r = range(v);

			Size gs, ls;
			size(v, gs);
			local_size(v, ls);
			VecSetType(result.implementation(), VECMPI);
			VecSetSizes(result.implementation(), ls.get(0), gs.get(0));

			writeLock(result);
			readLock(v);

			for(PetscInt i = r.begin(); i < r.end(); ++i) {
				PetscScalar value = get(v, i);
				VecSetValue(result.implementation(), i, -value, INSERT_VALUES);
			}

			writeUnlock(result);
			readUnlock(v);
			return true;
		}

		PetscScalar trace(const Matrix &mat);

		template<int FillType>
		bool apply_tensor_reduce(const PETScGenericMatrix<FillType> &mat, const Plus &, const int dim, Vector &result) {
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

		template<int FillType>
		bool apply_tensor_reduce(const PETScGenericMatrix<FillType> &mat, const Min &, const int dim, Vector &result) {
			PetscInt grows, gcols;
			MatGetSize(mat.implementation(), &grows, &gcols);

			if (dim == 1) {
				result.init({PETSC_DECIDE}, {grows});
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

		template<int FillType>
		bool apply_tensor_reduce(const PETScGenericMatrix<FillType> &mat, const Max &, const int dim, Vector &result) {
			PetscInt grows, gcols;
			MatGetSize(mat.implementation(), &grows, &gcols);

			if (dim == 1) {
				result.init({PETSC_DECIDE}, {grows});
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

		bool inverse(const Matrix &mat, Matrix &result);

		// bool gemm(const PetscScalar alpha, const Matrix &left, const Matrix &right,
		// 		  const bool transpose_left, const bool transpose_right, const PetscScalar beta, Matrix &result);


		bool waxpby(const PetscScalar a, const Vector &x, const PetscScalar &b, const Vector &y, Vector &result);
	};

	template<>
	class Backend<PetscScalar, PETSC> : public PETScBackend {
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
			info_.setName("petsc");
		}
	};
}

#endif //UTOPIA_UTOPIA_PETSCBACKEND_HPP

