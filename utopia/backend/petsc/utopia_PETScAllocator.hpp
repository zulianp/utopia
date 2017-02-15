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
		static std::unique_ptr<Mat, std::function<void(Mat*)>> claim(MPI_Comm comm, Size local, Size global) {
			assert(global.n_dims() == 2);
			assert(local.n_dims() == 2);

			//TODO ask pool for m
			Mat* m = new Mat;
			MatCreate(comm, m);

			return std::unique_ptr<Mat, std::function<void(Mat*)>>(m, mat_destructor);
		}

		static std::unique_ptr<Mat, std::function<void(Mat*)>> clone(const std::unique_ptr<Mat, std::function<void(Mat*)>>&) {
			assert(false && "TODO");
			return nullptr;
		}
	};

}

#endif // UTOPIA_PETSC_ALLOCATOR_HPP
#endif //WITH_PETSC
