#include "utopia_petsc_Eval_Blocks.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include <vector>

namespace utopia {

	void build_blocks(PetscMatrix &left, const Blocks<PetscMatrix> &blocks)
	{
	    std::vector<Mat> matrices;

	    MPI_Comm comm = PETSC_COMM_WORLD;
	    for(auto b_ptr : blocks.blocks()) {

	        if(b_ptr) {
	            comm = b_ptr->comm().get();
	            matrices.push_back(b_ptr->raw_type());
	        }  else {
	            matrices.push_back(nullptr);
	        }
	    }

        //FIXME the comm might be a sub-comm and needs managed memory
	    left.nest(comm, blocks.rows(), nullptr, blocks.cols(), nullptr, &matrices[0]);
	}

    void build_blocks(PetscVector &left, const Blocks<PetscVector> &blocks)
    {
        std::vector<Vec> vectors;

        MPI_Comm comm = PETSC_COMM_WORLD;
        for(auto b_ptr : blocks.blocks()) {

            if(b_ptr) {
                comm = b_ptr->comm().get();
                vectors.push_back(b_ptr->raw_type());
            }  else {
                vectors.push_back(nullptr);
            }
        }

        //FIXME the comm might be a sub-comm and needs managed memory
        left.nest(comm, blocks.size(), nullptr, &vectors[0]);
    }

}
