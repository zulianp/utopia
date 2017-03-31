
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H


#include "utopia_PETScError.hpp"
#include "utopia_Memory.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {

	template<>
	class Allocator<Mat, FillType::DENSE> {
	public:
		static void destructor(Mat* m) {
			MEMPOOL().put(m);
		};

		static MemoryPtr<Mat> claim(MPI_Comm comm, const Size& local, const Size& global) {
			MEMPOOL().setCommunicator(comm);
			Mat* m = MEMPOOL().getMat(local, global);

			return MemoryPtr<Mat>(m, destructor);
		}

		static MemoryPtr<Mat> clone(const MemoryPtr<Mat>& m) {
			Mat* new_m = MEMPOOL().getMat(*m);
			MatCopy(*m, *new_m, SAME_NONZERO_PATTERN);

			return MemoryPtr<Mat>(new_m, destructor);
		}
	};


	template<int FillType = FillType::DENSE>
	using PETScGenericMatrix = Memory<Mat, FillType>;

	using PETScMatrix = PETScGenericMatrix<>;

}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H
