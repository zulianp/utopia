
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

	static std::function<void(Vec*)> vec_destructor = [] (Vec* v) {
		//TODO return v to pool
		VecDestroy(v);
		delete v;
	};
	
	template<>
	class Allocator<Vec, 0> {
	public:
		static std::unique_ptr<Vec, std::function<void(Vec*)>> claim(MPI_Comm comm, Size local, Size global) {
			assert(global.n_dims() == 1);
			assert(local.n_dims() == 1);
	
			//TODO ask pool for v
			Vec* v = new Vec;
			VecCreateMPI(comm, local.get(0), global.get(0), v);
	
			return std::unique_ptr<Vec, std::function<void(Vec*)>>(v, vec_destructor);
		}
	
		static std::unique_ptr<Vec, std::function<void(Vec*)>> clone(const std::unique_ptr<Vec, std::function<void(Vec*)>>& v) {
			Vec* new_v = new Vec;
			VecDuplicate(*v, new_v);
			VecCopy(*v, *new_v);
	
			return std::unique_ptr<Vec, std::function<void(Vec*)>>(new_v, vec_destructor);
		}
	};


    // Newer PETScVector (one level wrapper)

    class PETScVector : public Memory<Vec, 0> {

    public:
        PETScVector(MPI_Comm comm = PETSC_COMM_WORLD) : super(comm) {
            Vec* v = new Vec;
            VecCreate(comm, v);
            mem_ = std::unique_ptr<Vec, std::function<void(Vec*)>>(v, vec_destructor);
        }

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
