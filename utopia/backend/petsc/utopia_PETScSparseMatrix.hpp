
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_PETScMatrix.hpp"

namespace utopia{

	// TODO - find a way to pass NNZ to the Allocator
	template<>
	class Allocator<Mat, FillType::SPARSE> {
	public:
		static void destructor(Mat* m) {
			//TODO return m to pool
			MatDestroy(m);
			delete m;
		};

		static MemoryPtr<Mat> claim(MPI_Comm comm, const Size& local, const Size& global) {
			//TODO ask pool for m
			Mat* m = new Mat;
			if (global.n_dims() >= 2 && local.n_dims() >= 2) {
				MatCreateAIJ(comm, local.get(0), local.get(1), global.get(0), global.get(1),
					1, PETSC_NULL, 1, PETSC_NULL, m);
			} else {
				MatCreate(comm, m);
			}

			return MemoryPtr<Mat>(m, destructor);
		}

		static MemoryPtr<Mat> clone(const MemoryPtr<Mat>& m) {
			Mat* new_m = new Mat;
			MatDuplicate(*m, MAT_COPY_VALUES, new_m);

			return MemoryPtr<Mat>(new_m, destructor);
		}
	};


	using PETScSparseMatrix = PETScGenericMatrix<FillType::SPARSE>;
}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
