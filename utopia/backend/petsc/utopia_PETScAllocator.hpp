#include "utopia_Memory.hpp"

#ifdef WITH_PETSC

#ifndef UTOPIA_PETSC_HPP
#define UTOPIA_PETSC_HPP

#include "petscmat.h"
#include "petscvec.h"


template<>
class Allocator<Mat, 0> {
public:
	static std::shared_ptr<Mat> claim(Size global, Size local) {
		assert(global.n_dims() == 2);
		assert(local.n_dims() == 2);

		//TODO ask pool for m
		Mat* m = new Mat;
		MatCreate(PETSC_COMM_WORLD, m);

		return std::shared_ptr<Mat>(m, [] (Mat* m){
			//TODO return m to pool
			MatDestroy(m);
			delete m;
		});
	}
};

template<>
class Allocator<Vec, 0> {
public:
	static std::shared_ptr<Vec> claim(Size global, Size local) {
		assert(global.n_dims() == 1);
		assert(local.n_dims() == 1);

		//TODO ask pool for v
		Vec* v = new Vec;
		VecCreateMPI(PETSC_COMM_WORLD, local.get(0), global.get(0), v);

		return std::shared_ptr<Vec>(v, [] (Vec* v){
			//TODO return v to pool
			VecDestroy(v);
			delete v;
		});
	}
};


#endif // UTOPIA_PETSC_HPP
#endif //WITH_PETSC
