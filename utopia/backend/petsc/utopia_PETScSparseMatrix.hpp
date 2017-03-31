
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_PETScMatrix.hpp"

namespace utopia{

	// TODO - find a way to pass NNZ to the Allocator
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
			MatConvert(*m, MATAIJ, MAT_INITIAL_MATRIX, new_m);

			return MemoryPtr<Mat>(new_m, destructor);
		}
	};


	using PETScSparseMatrix = PETScGenericMatrix<FillType::SPARSE>;
}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
