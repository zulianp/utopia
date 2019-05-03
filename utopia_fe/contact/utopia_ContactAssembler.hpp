#ifndef UTOPIA_CONTACT_ASSEMBLER_HPP
#define UTOPIA_CONTACT_ASSEMBLER_HPP

#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_Contact.hpp"

#include <cassert>

namespace utopia {

   class ContactAssembler {
   public:
    bool assemble(
        libMesh::MeshBase &mesh,
        libMesh::DofMap &dof_map,
        const ContactParams &params);
    };

}

#endif //UTOPIA_CONTACT_ASSEMBLER_HPP