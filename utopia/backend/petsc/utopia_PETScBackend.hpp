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
		bool convert(PETScMatrix &wrapper, Mat mat);
		bool convert(PETScSparseMatrix &wrapper, Mat mat);
		bool wrap(Mat mat, PETScSparseMatrix &wrapper);
		bool wrap(Vec vec, PETScVector &wrapper);
		
		Range range(const PETScVector &v);
		Range rowRange(const PETScMatrix &m);
		Range colRange(const PETScMatrix &m);
		
		bool size(const PETScMatrix &m, Size &size) const;
		bool size(const PETScVector &m, Size &size) const;
		bool local_size(const PETScVector &m, Size &size) const;
		bool local_size(const PETScMatrix &m, Size &size) const;

		void resize(const Size &size, PETScMatrix &mat);
		void resize(const Size &local_s, const Size &global_s, PETScMatrix &mat);
		void resize(const Size &size, PETScVector &mat);
		void resize(const Size &s_local, const Size &s_global, PETScVector &vec);
		
		void assignFromRange(PETScMatrix &left, const PETScMatrix &right, const Range &globalRowRange,
							 const Range &globalColRange);
		
		// read matrix
		bool read(const std::string &path, PETScMatrix &Mat_A);
		// write matrix
		bool write(const std::string &path, const PETScMatrix &Mat_A);
		
		// write vector
		bool write(const std::string &path, const PETScVector &Vec_A);
		
		// monitor for cyrill 
		bool monitor(const long & it,  PETScMatrix &Mat_A); 
		bool monitor(const long & it,  PETScVector &Vec_A); 


		PetscScalar get_global_nnz(PETScMatrix &Mat_A); 
		PetscScalar get_local_nnz(PETScMatrix &Mat_A); 

		// call to MatZeroRows
		bool set_zero_rows(PETScMatrix &Mat_A, const std::vector<int> &index); 
		bool apply_BC_to_system(PETScMatrix & A, PETScVector& x, PETScVector& rhs, const std::vector<int> &index); 

		// read vector
		bool read(const std::string &path, PETScVector &Vec_A);
		
		void assignFromRange(PETScVector &left, const PETScVector &right, const Range &globalRowRange,
							 const Range & /*globalColRange */);
		
		bool scal(const PetscScalar /*scaleFactor*/, const PETScMatrix & /*m*/, PETScMatrix &/*result*/);
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
				std::cerr << "[Error] axpy" << std::endl;
			}
			
			return aux_zaxpy(scaleFactor, left, result);
		}
		
		bool aux_zaxpy(const PetscScalar scaleFactor, const PETScVector &left,
					   PETScVector &result);
		
		bool aux_zaxpy(const PetscScalar scaleFactor, const PETScMatrix &left,
					   PETScMatrix &result);
		
		bool gemm(const PetscScalar /*alpha*/, const PETScMatrix &/*left*/, const PETScMatrix &/*right*/, const PetscScalar /* beta*/,
				  PETScMatrix & /*result*/);
		
		bool gemm(const PetscScalar /*alpha*/, const PETScMatrix &/*left*/, const PETScMatrix &/*right*/,
				  bool /*transpose_left*/, bool /*transpose_right*/, const PetscScalar /* beta*/, PETScMatrix &/*result*/);
		
		bool gemm(const PetscScalar /*alpha */, const PETScVector &/*left*/, const PETScMatrix &/*right*/,
				  bool /*transpose_left*/, bool /*transpose_right*/, const PetscScalar /*beta*/, PETScMatrix &/*result*/);
		
		
		template<class Left, class Right, class Result>
		bool gem(const double /*scaleFactor*/, const Left &/*left*/, const Right &/*right*/, bool /*transpose_left*/,
				 bool /*transpose_right*/, const double /*beta*/, Result &/*result*/);

  		bool gemm(const double scaleFactor, const PETScMatrix &left, const Vector &right, bool transpose_left,  bool transpose_right, const double beta, Vector &result);
  

		bool transpose(const PETScMatrix &mat, PETScMatrix &result);
		
		bool scal(const PetscScalar /*scale_factor*/, const PETScVector &/*v*/, PETScVector &/*result */);
		
		void build(PETScMatrix &m, const Size &size, const Identity &);

		// inline void init(PETScMatrix &mat)
		// {

		// }
		
		template<class Tensor>
		void build(Tensor &t, const Size &s, const Resize &)
		{
			//FIXME switch with more efficient version without zeros
			build(t, s, Zeros());
			// resize(s, t);
		}
		
		void build(PETScSparseMatrix &m, const Size &size, const Identity &);
		
		void build(PETScMatrix &m, const Size &size, const LocalIdentity &);
		
		void build(PETScSparseMatrix &m, const Size &size, const LocalIdentity &);
		
		void build(PETScSparseMatrix &m, const Size &size, const NNZ<PetscInt> &nnz);
		
		void build(PETScSparseMatrix &m, const Size &size, const LocalNNZ<PetscInt> &nnz);
		void build(PETScSparseMatrix &m, const Size &size, const LocalRowNNZ<PetscInt> &nnz);
		/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
		void build(PETScMatrix  &m, const Size &size, const LocalNNZ<PetscInt> & /*nnz */);
		
		/// Obviously there is no sparse support for dense matrices. Nevertheless, compatibility requires it.
		void build(PETScMatrix  &m, const Size &size, const NNZ<PetscInt> &/*nnz*/);
		
		void build(PETScMatrix &m, const Size &size, const Zeros &);
		
		void build(PETScVector &v, const Size &size, const Zeros &);
		
		void build(PETScMatrix &m, const Size &size, const LocalZeros &);
		void build(PETScVector &v, const Size &size, const LocalZeros &);
		
		void build(PETScMatrix &m, const Size &size, const Values<PetscScalar> &values);
		
		void build(PETScVector &v, const Size &local_size, const Size &&global_size, const Values<PetscScalar> &values);
		
		void build(PETScVector &v, const Size &size, const Values<PetscScalar> &values);
		
		void build(PETScMatrix &m, const Size &size, const LocalValues<PetscScalar> &values);
		
		void build(PETScVector &v, const Size &size, const LocalValues<PetscScalar> &values);
		
		void set(PETScVector &v, const PetscInt index, PetscScalar value);
		void add(PETScVector &v, const PetscInt index, PetscScalar value);
		
		void set(PETScVector &v, const std::vector<PetscInt> indices, const std::vector<PetscScalar> values);
		void set(PETScMatrix &v, const PetscInt row, const PetscInt col, PetscScalar value);
		void set(PETScSparseMatrix &v, const PetscInt row, const PetscInt col, PetscScalar value);
		
		void add(PETScMatrix &v, const PetscInt row, const PetscInt col, PetscScalar value);
		
		template<class Tensor>
		void readLock(const Tensor &) {}
		
		template<class Tensor>
		void readUnlock(const Tensor &) {}
		
		void writeLock(PETScVector &vec);
		void writeUnlock(PETScVector &vec);
		void writeLock(const PETScMatrix &mat);
		void writeUnlock(const PETScMatrix &mat);
		
		void writeLock(const PETScSparseMatrix &mat);
		void writeUnlock(const PETScSparseMatrix &mat);
		
		
		void set(PETScMatrix &v, const std::vector<PetscInt> rows, const std::vector<PetscInt> cols, const std::vector<PetscScalar> values);
		PetscScalar get(const PETScVector &v, const PetscInt index);
		PetscScalar get(const PETScMatrix &v, const PetscInt row, const PetscInt col);
		bool apply(const PETScMatrix &left, const PETScVector &right, const Multiplies &, PETScVector &result);
		bool apply(const PETScMatrix &left, const PETScMatrix &right, const Multiplies &, PETScMatrix &result);
		bool apply(const PETScVector &left, const PETScVector &right, const EMultiplies &, PETScVector &result);
		bool apply(const PETScVector &left, const PETScVector &right, const Divides &, PETScVector &result);
		
		bool mat_mult_add(const PETScMatrix &m, const PETScVector &right, const PETScVector &left, PETScVector &result);
		bool mat_multT_add(const PETScMatrix &m, const PETScVector &right, const PETScVector &left, PETScVector &result);
		
		// reciprocal
		bool apply(const PETScVector &vec, const Reciprocal<double> &reciprocal, PETScVector &result);
		PetscScalar norm2(const PETScVector &v);
		PetscScalar norm2(const Matrix &m);
		PetscScalar norm1(const PETScVector &v);
		PetscScalar norm_infty(const PETScVector &v);
		PetscScalar reduce(const PETScVector &vec, const Plus &);
		PetscScalar reduce(const PETScMatrix &mat, const Plus &op);
		PetscScalar reduce(const PETScVector &, const Min &);
		PetscScalar reduce(const PETScMatrix &, const Min &);
		PetscScalar reduce(const PETScVector &, const Max &);
		PetscScalar reduce(const PETScMatrix &, const Max &);
		
		// get diagonal of matrix as vector
		bool diag(PETScVector &vec, const PETScMatrix &mat);
		
		bool diag(PETScSparseMatrix &mat, const PETScVector &vec);
		bool diag(PETScMatrix &mat, const PETScVector &vec);
		bool diag(PETScMatrix &out, const PETScMatrix &in);

        bool mat_diag_shift(PETScMatrix &left, const PetscScalar diag_factor);
        
		bool compare(const Vector &left, const Vector &right, const ApproxEqual &comp);
		
		bool compare(const Matrix &left, const Matrix &right, const ApproxEqual &comp);
		
		
		bool apply(const PetscScalar factor, const Vector &vec, const Multiplies &, Vector &result);
		bool apply(const PetscScalar factor, const Matrix &mat, const Multiplies &, Matrix &result);

		bool apply(const Vector &vec, const Abs &, Vector &result);
		// bool apply(const Matrix &mat, const Abs &, Matrix &result);
		
		PetscScalar dot(const Vector &left, const Vector &right) const;
		
		
		bool mul(const PETScMatrix &left, const PETScMatrix &right, PETScMatrix &result);
		
		bool outer(const Vector &left, const Vector &right, Matrix &result);
		
		void assignTransposed(PETScMatrix &left, const PETScMatrix &right);
		
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

		bool triple_product_PtAP(const PETScMatrix &, const PETScMatrix &, PETScMatrix &); 
		bool triple_product(const PETScMatrix &, const PETScMatrix &, const PETScMatrix &, PETScMatrix &); 

		bool is_nan_or_inf(const PETScVector &v); 
		bool is_nan_or_inf(const PETScMatrix &m); 


		template<class Tensor>
		void readAndWriteLock(Tensor &t) {
			writeLock(t);
		}
		
		template<class Tensor>
		void readAndWriteUnlock(Tensor &t){
			writeUnlock(t);
		}
		
		template<class LeftTensor, class RightTensor, class ResultTensor>
		inline bool apply(LeftTensor &&left, RightTensor &&right, const Plus &, ResultTensor &result) {
			return zaxpy(1.0, std::forward<LeftTensor>(left), std::forward<RightTensor>(right), result);
		}
		
		template<class LeftTensor, class RightTensor, class ResultTensor>
		inline bool apply(LeftTensor &&left, RightTensor &&right, const Minus &, ResultTensor &result) {
			return zaxpy(-1.0, std::forward<RightTensor>(right), std::forward<LeftTensor>(left), result);
		}

		bool diag_scale_right(const Matrix &m, const Vector &diag, Matrix &result);
		bool diag_scale_left(const Vector &diag, const Matrix &m,  Matrix &result);
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

		bool apply_tensor_reduce(const Matrix &mat, const Plus &, const int dim, Vector &result);
		bool apply_tensor_reduce(const Matrix &mat, const Min &, const int dim, Vector &result);
		bool apply_tensor_reduce(const Matrix &mat, const Max &, const int dim, Vector &result);
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

