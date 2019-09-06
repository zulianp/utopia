#ifndef UTOPIA_ADAPTIVITY_HPP
#define UTOPIA_ADAPTIVITY_HPP

#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "libmesh/dof_map.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/elem.h"

namespace utopia {

    class Adaptivity {
    public:
        void constraint_matrix(const LibMeshFunctionSpace &V, USparseMatrix &M, USparseMatrix &S);
        void constraint_matrix(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num, USparseMatrix &M, USparseMatrix &S);


    private:
        libMesh::DofConstraints dof_constraints_;

        // void assemble_constraint(const LibMeshFunctionSpace &V);
        void assemble_constraint(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num);

    public:
        static void compute_constraints(libMesh::DofConstraints &constraints,
                                        const libMesh::DofMap &dof_map,
                                        const unsigned int variable_number,
                                        const libMesh::Elem * elem,
                                        const unsigned mesh_dim
                                        );

        static void compute_all_constraints(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map,
            libMesh::DofConstraints &constraints);

        template<class ElementMatrix>
        static void construct_constraint_matrix(const libMesh::Elem *elem,
                                                const libMesh::DofMap &dof_map,
                                                const libMesh::DofConstraints &dof_constraints,
                                                ElementMatrix &mat,
                                                std::vector<libMesh::dof_id_type> &dof_indices,
                                                const bool called_recursively = false)
        {
            const std::size_t n_var = dof_map.n_variables();
            dof_map.dof_indices(elem, dof_indices);

            typedef std::set<libMesh::dof_id_type> RCSet;

            RCSet dof_set;

            for(const auto &dof : dof_indices){
                if(dof_map.is_constrained_dof(dof)) {
                    auto pos = dof_constraints.find(dof);
                    if(pos == dof_constraints.end()) continue;

                    const auto &constraint_row = pos->second;

                    for(const auto &item : constraint_row) {
                        dof_set.insert(item.first);
                    }
                }
            }

            for(const auto & dof : dof_indices) {
                dof_set.erase(dof);
            }

            if (!dof_set.empty() ||  // case 1: constrained in terms of other DOFs
                !called_recursively) // case 2: constrained in terms of our own DOFs
            {
                const unsigned int old_size =
                static_cast<unsigned int>(dof_indices.size());

                // Add new dependency dofs to the end of the current dof set
                dof_indices.insert(dof_indices.end(),
                                   dof_set.begin(), dof_set.end());

                mat.resize(old_size, static_cast<unsigned int>(dof_indices.size()));
                mat.zero();

                for(unsigned int i=0; i != old_size; i++) {
                    if (dof_map.is_constrained_dof(dof_indices[i])) {
                        auto pos = dof_constraints.find(dof_indices[i]);

                        if(pos == dof_constraints.end()) {
                            mat(i, i) = 1.;
                            continue;
                        }

                        const auto & constraint_row = pos->second;

                        for(const auto & item : constraint_row) {
                            const auto n_elem_dofs = static_cast<unsigned int>(dof_indices.size());

                            for (unsigned int j=0; j != n_elem_dofs; j++) {
                                if (dof_indices[j] == item.first){
                                    mat(i, j) = item.second;
                                }
                            }
                        }

                    } else {
                        mat(i,i) = 1.;
                    }
                }
            }
        }


        template<class ElementMatrix>
        static void constrain_matrix(const libMesh::Elem *elem,
                                     const libMesh::DofMap &dof_map,
                                     const libMesh::DofConstraints dof_constraints,
                                     ElementMatrix &mat,
                                     std::vector<libMesh::dof_id_type> &dof_indices)
        {
            using uint = unsigned int;

            ElementMatrix C;
            
            construct_constraint_matrix(
                                        elem,
                                        dof_map,
                                        dof_constraints,
                                        C,
                                        dof_indices,
                                        false //FIXME
                                        );

            mat.left_multiply_transpose(C);
            mat.right_multiply(C);

            const uint ndofs = dof_indices.size();

            for(uint i = 0; i < ndofs; ++i) {
                if(dof_map.is_constrained_dof(dof_indices[i])) {
                    auto it = dof_constraints.find(dof_indices[i]);
                    if(it == dof_constraints.end()) continue;

                    for(uint j = 0; j < ndofs; j++) {
                        mat(i, j) = 0;
                    }

                    mat(i, i) = 1;


                    const auto &c_row = it->second;

                    for(const auto &item : c_row) {

                        for(uint j = 0; j < ndofs; j++) {
                            if(dof_indices[j] == item.first) {
                                mat(i, j) = -item.second;
                            }
                        }
                    }
                }
            }
        }

        template<class ElementVector>
        static void constrain_vector(const libMesh::Elem *elem,
                                     const libMesh::DofMap &dof_map,
                                     const libMesh::DofConstraints dof_constraints,
                                     ElementVector &vec,
                                     std::vector<libMesh::dof_id_type> &dof_indices)
        {
            using uint = unsigned int;

            libMesh::DenseMatrix<double> C;

            construct_constraint_matrix(
                                        elem,
                                        dof_map,
                                        dof_constraints,
                                        C,
                                        dof_indices,
                                        false //FIXME
                                        );

            libMesh::DenseVector<double> old(vec);
            C.vector_mult_transpose(vec, old);
        }

        template<class ElementMatrix, class ElementVector>
        static void constrain_matrix_and_vector(
            const libMesh::Elem *elem,
            const libMesh::DofMap &dof_map,
            const libMesh::DofConstraints dof_constraints,
            ElementMatrix &mat,
            ElementVector &vec,
            std::vector<libMesh::dof_id_type> &dof_indices)
        {
            using uint = unsigned int;

            libMesh::DenseMatrix<double> C;
            construct_constraint_matrix(
                                        elem,
                                        dof_map,
                                        dof_constraints,
                                        C,
                                        dof_indices,
                                        false //FIXME
                                        );

            libMesh::DenseVector<double> old(vec);
            C.vector_mult_transpose(vec, old);

            mat.left_multiply_transpose(C);
            mat.right_multiply(C);

            const uint ndofs = dof_indices.size();

            for(uint i = 0; i < ndofs; ++i) {
                if(dof_map.is_constrained_dof(dof_indices[i])) {
                    auto it = dof_constraints.find(dof_indices[i]);
                    if(it == dof_constraints.end()) continue;

                    for(uint j = 0; j < ndofs; j++) {
                        mat(i, j) = 0;
                    }

                    mat(i, i) = 1;

                    const auto &c_row = it->second;

                    for(const auto &item : c_row) {

                        for(uint j = 0; j < ndofs; j++) {
                            if(dof_indices[j] == item.first) {
                                mat(i, j) = -item.second;
                            }
                        }
                    }
                }
            }
        }
    };
}

#endif
