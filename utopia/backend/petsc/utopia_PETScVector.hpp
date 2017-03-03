
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

	template<>
	class Allocator<Vec, 0> {
	public:
		static void destructor(Vec* v) {
			//TODO return v to pool
			VecDestroy(v);
			delete v;
		};

		static MemoryPtr<Vec> claim(MPI_Comm comm, const Size& local, const Size& global) {
			//TODO ask pool for v
			Vec* v = new Vec;
			if (global.n_dims() >= 1 && local.n_dims() >= 1) {
				VecCreateMPI(comm, local.get(0), global.get(0), v);
			} else {
				VecCreate(comm, v);
			}

			return MemoryPtr<Vec>(v, destructor);
		}

		static MemoryPtr<Vec> clone(const MemoryPtr<Vec>& v) {
			Vec* new_v = new Vec;
			VecDuplicate(*v, new_v);
			VecCopy(*v, *new_v);

			return MemoryPtr<Vec>(new_v, destructor);
		}
	};

	typedef Memory<Vec, 0> PETScVector;

}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
