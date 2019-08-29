#ifndef UTOPIA_ADAPTIVITY_HPP
#define UTOPIA_ADAPTIVITY_HPP

#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "libmesh/dof_map.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

namespace utopia {

    class Adaptivity {
    public:
        void constraint_matrix(const LibMeshFunctionSpace &V, USparseMatrix &M);
        void constraint_matrix(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num, USparseMatrix &M);

    private:
        libMesh::DofConstraints dof_constraints_;

        // void assemble_constraint(const LibMeshFunctionSpace &V);
        void assemble_constraint(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num);
        
        void compute_constraints(
            libMesh::DofConstraints &constraints, 
            const libMesh::DofMap &dof_map,
            const unsigned int variable_number, 
            const libMesh::Elem * elem,
            const unsigned mesh_dim
        );
    };
}

#endif