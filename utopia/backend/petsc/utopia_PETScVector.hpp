
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
			MEMPOOL().put(v);
		};

		static MemoryPtr<Vec> claim(MPI_Comm comm, const Size& local, const Size& global) {
			MEMPOOL().setCommunicator(comm);
			Vec* v = MEMPOOL().getVec(local.n_dims() > 1 ? Size({local.get(0)}) : local,
				global.n_dims() > 1 ? Size({global.get(0)}) : global);

			return MemoryPtr<Vec>(v, destructor);
		}

		static MemoryPtr<Vec> clone(const MemoryPtr<Vec>& v) {
			Vec* new_v = MEMPOOL().getVec(*v);
			VecCopy(*v, *new_v);

			return MemoryPtr<Vec>(new_v, destructor);
		}
	};

	typedef Memory<Vec, 0> PETScVector;

}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
