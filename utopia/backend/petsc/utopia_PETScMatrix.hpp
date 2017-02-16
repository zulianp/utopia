
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H


#include "utopia_PETScError.hpp"
#include "utopia_PETScVector.hpp"
#include "utopia_Base.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {

	static void mat_destructor(Mat* m){
		//TODO return m to pool
		MatDestroy(m);
		delete m;
	};

	template<>
	class Allocator<Mat, 0> {
	public:
		static MemoryPtr<Mat> claim(MPI_Comm comm, const Size& local, const Size& global) {
			//TODO ask pool for m
			Mat* m = new Mat;
			if (global.n_dims() == 2 && local.n_dims() == 2) {
				assert(false);
			} else {
				MatCreate(comm, m);
			}

			return MemoryPtr<Mat>(m, mat_destructor);
		}

		static MemoryPtr<Mat> clone(const MemoryPtr<Mat>& m) {
			Mat* new_m = new Mat;
			MatDuplicate(*m, MAT_COPY_VALUES, new_m);

			return MemoryPtr<Mat>(new_m, mat_destructor);
		}
	};


	class PETScMatrix : public Memory<Mat, 0> {
	public:
		using Memory<Mat, 0>::Memory;

		void describe() const {
			MatView(Memory<Mat, 0>::implementation(), PETSC_VIEWER_STDOUT_WORLD);
		}
	
		// void save(const std::string &path) {
		//     PetscViewer viewer;
		//
		//     PetscViewerCreate(_wrapper->communicator(), &viewer);
		//     PetscViewerSetType(viewer, PETSCVIEWERASCII);
		//     PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
		//     PetscViewerFileSetName(viewer, path.c_str());
		//     PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
		//     MatView(_wrapper->implementation(), viewer);
		//     PetscViewerDestroy(&viewer);
		// }
		//
		// void setName(const std::string &name) {
		//     PetscObjectSetName((PetscObject) _wrapper->implementation(), name.c_str());
		// }

	};
}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H
