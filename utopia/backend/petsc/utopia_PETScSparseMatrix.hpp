
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_PETScMatrix.hpp"

namespace utopia{

	template<>
	inline void check<Mat, FillType::SPARSE>(Mat m) {
		// std::cout << "called spec sparse: ";
		MatType type;
		if (MatGetType(m, &type) == 0) {
			if (type && strcmp(MATSEQAIJ, type) && strcmp(MATMPIAIJ, type)) {
				std::cout << "[Error] Placing a dense matrix in a sparse one!" << '\n';
				std::cout << "        Matrix is of type " << type << '\n';
				assert(false && "Cannot wrap a dense matrix in a sparse wrapper");
			}
		}
		// std::cout << type << '\n';
	}

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
			Mat* new_m = MEMPOOL().getSparseMat(*m);
			// MatSeqAIJSetPreallocation(*new_m, m_info.nz_used, NULL); // this works in tests
			MatSetOption(*new_m, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); // but this is way faster
			MatCopy(*m, *new_m, DIFFERENT_NONZERO_PATTERN);

			return MemoryPtr<Mat>(new_m, destructor);
		}
	};


	using PETScSparseMatrix = PETScGenericMatrix<FillType::SPARSE>;
}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
