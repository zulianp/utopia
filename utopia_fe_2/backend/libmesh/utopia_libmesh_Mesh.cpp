#include "utopia_libmesh_Mesh.hpp"

// utopia
#include "utopia_make_unique.hpp"

// utopia fe

// libmesh
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel_mesh.h"

namespace utopia {

    Mesh<libMesh::MeshBase>::Mesh(const Communicator &comm)
        : comm_(std::make_shared<libMesh::Parallel::Communicator>(comm.raw_comm())),
          impl_(utopia::make_unique<libMesh::DistributedMesh>(*comm_)) {}

    Mesh<libMesh::MeshBase>::~Mesh() {}

    void Mesh<libMesh::MeshBase>::read(Input &) {}
}  // namespace utopia