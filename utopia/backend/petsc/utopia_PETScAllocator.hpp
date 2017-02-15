#ifdef WITH_PETSC

#ifndef UTOPIA_PETSC_ALLOCATOR_HPP
#define UTOPIA_PETSC_ALLOCATOR_HPP

#include "utopia_Memory.hpp"

#include "petscmat.h"
#include "petscvec.h"

namespace utopia {

	static std::function<void(Mat*)> mat_destructor = [] (Mat* m){
		//TODO return m to pool
		MatDestroy(m);
		delete m;
	};

	template<>
	class Allocator<Mat, 0> {
	public:
		static MemoryPtr<Mat> claim(MPI_Comm comm, const Size& local, const Size& global) {
			assert(global.n_dims() == 2);
			assert(local.n_dims() == 2);

			//TODO ask pool for m
			Mat* m = new Mat;
			MatCreate(comm, m);

			return MemoryPtr<Mat>(m, mat_destructor);
		}

		static MemoryPtr<Mat> clone(const MemoryPtr<Mat>&) {
			assert(false && "TODO");
			return nullptr;
		}
	};

}

#endif // UTOPIA_PETSC_ALLOCATOR_HPP
#endif //WITH_PETSC
