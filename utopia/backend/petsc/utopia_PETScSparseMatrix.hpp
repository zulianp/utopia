
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_PETScMatrix.hpp"

namespace utopia{

	// NNZ has to be set in the backend, after allocation of the Matrix
	template<>
	class Allocator<Mat, FillType::SPARSE> {
	public:
		static void destructor(Mat* m) {
			MEMPOOL().putSparse(m);
		};

		static MemoryPtr<Mat> claim(MPI_Comm comm, const Size& local, const Size& global) {
			MEMPOOL().setCommunicator(comm);
			Mat* m = MEMPOOL().getSparseMat(local, global);

			return MemoryPtr<Mat>(m, destructor);
		}

		static MemoryPtr<Mat> clone(const MemoryPtr<Mat>& m) {
			// FIXME - cannot find a way to reuse sparse matrices with MatCopy
			Mat* new_m = new Mat; //MEMPOOL().getSparseMat(*m);
			MatDuplicate(*m, MAT_COPY_VALUES, new_m);

			return MemoryPtr<Mat>(new_m, destructor);
		}
	};


	using PETScSparseMatrix = PETScGenericMatrix<FillType::SPARSE>;
}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
