#ifndef MIXED_PAR_MORTAR_ASSEMBLER_HPP
#define MIXED_PAR_MORTAR_ASSEMBLER_HPP

#include <memory>
#include <vector>
#include "utopia_fe.hpp"
#include "utopia_LibMeshBackend.hpp"
#include <libmesh/sparse_matrix.h>
#include "Box.hpp"
#include "libmesh/serial_mesh.h"

namespace utopia {

	class MixedParMortarAssembler {
	public:

		MixedParMortarAssembler(
           libMesh::Parallel::Communicator &libmesh_comm,
           // const MPI_Comm comm,
			const std::shared_ptr<LibMeshFESpaceBase> &master,
			const std::shared_ptr<LibMeshFESpaceBase> &slave); //constructor
        
        void set_use_biorthogonal_multipliers(const bool use_biorth)
        {
            use_biorth_ = use_biorth;
        }
        
        bool Assemble(DSMatrixd &B);
        bool Transfer(DSMatrixd &B,DSMatrixd &T);


    
    private:
        libMesh::Parallel::Communicator &libmesh_comm_;
       // MPI_Comm comm_;
		std::shared_ptr<LibMeshFESpaceBase> master_;
		std::shared_ptr<LibMeshFESpaceBase> slave_;
        std::shared_ptr<LibMeshFESpaceBase> space_;
        bool use_biorth_;
		//std::vector< std::shared_ptr<MortarIntegrator> > integrators_;
	};
    


}

#endif //MFEM_L2P_PAR_MORTAR_ASSEMBLER_HPP
