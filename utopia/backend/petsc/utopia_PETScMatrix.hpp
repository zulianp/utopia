
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H


#include "utopia_PETScError.hpp"
#include "utopia_PETScVector.hpp"
#include "utopia_Base.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {

	template<>
	class Allocator<Mat, 0> {
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
				MatCreateDense(comm, local.get(0), local.get(1), global.get(0), global.get(1), NULL, m);
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

	typedef Memory<Mat, 0> PETScMatrix;

}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H
