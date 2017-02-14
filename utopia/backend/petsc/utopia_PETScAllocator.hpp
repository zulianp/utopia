#include "utopia_Memory.hpp"

#ifdef WITH_PETSC

#ifndef UTOPIA_PETSC_HPP
#define UTOPIA_PETSC_HPP

#include "petscmat.h"
#include "petscvec.h"

namespace utopia {

	template<>
	class Allocator<Mat, 0> {
	public:
		static std::unique_ptr<Mat> claim(Size global, Size local) {
			assert(global.n_dims() == 2);
			assert(local.n_dims() == 2);

			//TODO ask pool for m
			Mat* m = new Mat;
			MatCreate(PETSC_COMM_WORLD, m);

			return std::unique_ptr<Mat>(m, [] (Mat* m){
				//TODO return m to pool
				MatDestroy(m);
				delete m;
			});
		}

		static std::unique_ptr<Mat> clone(std::unique_ptr<Mat>&) {
			assert(false && "TODO");
			return nullptr;
		}
	};

	template<>
	class Allocator<Vec, 0> {
	public:
		static std::unique_ptr<Vec> claim(Size global, Size local) {
			assert(global.n_dims() == 1);
			assert(local.n_dims() == 1);

			//TODO ask pool for v
			Vec* v = new Vec;
			VecCreateMPI(PETSC_COMM_WORLD, local.get(0), global.get(0), v);

			return std::unique_ptr<Vec>(v, destructor_);
		}

		static std::unique_ptr<Vec> clone(std::unique_ptr<Vec>& v) {
			Vec* new_v = new Vec;
			VecDuplicate(*v, new_v);
			VecCopy(*v, *new_v);

			return std::unique_ptr<Vec>(new_v, destructor_);
		}

	private:
		std::function destructor_ = [] (Vec* v){
			//TODO return v to pool
			VecDestroy(v);
			delete v;
		};
	};

}

#endif // UTOPIA_PETSC_HPP
#endif //WITH_PETSC
