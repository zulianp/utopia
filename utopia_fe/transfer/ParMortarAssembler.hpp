// #ifndef PAR_MORTAR_ASSEMBLER_HPP
// #define PAR_MORTAR_ASSEMBLER_HPP

// #include "utopia_fe_core.hpp"
// #include "utopia_LibMeshBackend.hpp"
// #include "Box.hpp"
// #include "MortarAssembler.hpp"
// #include "moonolith_communicator.hpp"

// #include <libmesh/sparse_matrix.h>
// #include "libmesh/serial_mesh.h"

// #include <memory>
// #include <vector>

// namespace utopia {

// 	class ParMortarAssembler {
// 	public:

// 		ParMortarAssembler(
//             libMesh::Parallel::Communicator &libmesh_comm,
// 			const std::shared_ptr<LibMeshFESpaceBase> &master_slave); //constructor
        
        
//         bool Assemble(DSMatrixd &B);
//         bool SurfaceAssemble(DSMatrixd &B, const libMesh::Real search_radius, const int tag_1, const int tag_2);
//         bool Transfer(DSMatrixd &B,DSMatrixd &T);

    
//     private:
//         libMesh::Parallel::Communicator &libmesh_comm_;
// 		std::shared_ptr<LibMeshFESpaceBase> master_slave_;
//         std::shared_ptr<LibMeshFESpaceBase> space_;
// 		//std::vector< std::shared_ptr<MortarIntegrator> > integrators_;
// 	};
 
// }

// #endif //MFEM_L2P_PAR_MORTAR_ASSEMBLER_HPP
