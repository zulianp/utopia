
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H


#include "utopia_PETScError.hpp"
#include "utopia_Memory.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {

	template<>
	inline void check<Mat, FillType::DENSE>(Mat m) {
		// std::cout << "called spec dense: ";
		MatType type;
		if (MatGetType(m, &type) == 0) {
			if (type && !strcmp(MATSEQAIJ, type)) {
				std::cout << "[Error] Placing a sparse matrix in a dense one!" << '\n';
				assert(false && "Cannot wrap a sparse matrix in a dense wrapper");
			}
		}
		// std::cout << type << '\n';
	}

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
			MatType t;
			MatGetType(*m, &t);

			if (strcmp(MATSEQAIJ, t) == 0) {
				std::cout << "[Warning] Cloning sparse matrix like a dense one!" << '\n';
				assert(false);
			}
			Mat* new_m = MEMPOOL().getMat(*m);
			MatCopy(*m, *new_m, DIFFERENT_NONZERO_PATTERN);

			return MemoryPtr<Mat>(new_m, destructor);
		}
	};


	template<int FillType = FillType::DENSE>
	using PETScGenericMatrix = Memory<Mat, FillType>;

	using PETScMatrix = PETScGenericMatrix<>;

}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H
