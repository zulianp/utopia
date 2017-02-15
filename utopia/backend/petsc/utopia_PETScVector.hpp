
#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_PETScError.hpp"
#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Memory.hpp"

#include "petscmat.h"
#include "petscvec.h"

#include <memory>

namespace utopia {

	static void vec_destructor(Vec* v) {
		//TODO return v to pool
		VecDestroy(v);
		delete v;
	};

	template<>
	class Allocator<Vec, 0> {
	public:
		static MemoryPtr<Vec> claim(MPI_Comm comm, const Size& local, const Size& global) {
			//TODO ask pool for v
			Vec* v = new Vec;
			if (global.n_dims() == 1 && local.n_dims() == 1) {
				VecCreateMPI(comm, local.get(0), global.get(0), v);
			} else {
				VecCreate(comm, v);
			}

			return MemoryPtr<Vec>(v, vec_destructor);
		}

		static MemoryPtr<Vec> clone(const MemoryPtr<Vec>& v) {
			Vec* new_v = new Vec;
			VecDuplicate(*v, new_v);
			VecCopy(*v, *new_v);

			return MemoryPtr<Vec>(new_v, vec_destructor);
		}
	};


	// Newer PETScVector (one level wrapper)

	class PETScVector : public Memory<Vec, 0> {

	public:
		PETScVector(MPI_Comm comm = PETSC_COMM_WORLD) : super(comm) { }

		PETScVector(const PETScVector &other) : super(other) { }

		PETScVector &operator=(const PETScVector &other) {
			super::operator= (other);
			return *this;
		}

		~PETScVector() {
			destroy();
		}

		bool initialize(const PetscInt nLocal, const PetscInt nGlobal) {
			// FIXME Error check?
			super::initialize({nLocal}, {nGlobal});
			return true;
		}

		void resize(const PetscInt nLocal, const PetscInt nGlobal = PETSC_DETERMINE) {
			// FIXME Error check?
			super::resize({nLocal}, {nGlobal});
		}

		void describe() const {
			VecView(super::implementation(), PETSC_VIEWER_STDOUT_WORLD);
		}

	private:
		typedef Memory<Vec, 0> super;
	};
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
