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
        
        void constraint_matrix(const libMesh::MeshBase &mesh, 
                               const libMesh::DofMap &dof_map, 
                               USparseMatrix &M, USparseMatrix &S);


    private:
        libMesh::DofConstraints dof_constraints_;

        // void assemble_constraint(const LibMeshFunctionSpace &V);
        void assemble_constraint(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map);

    public:
        static void compute_constraints(libMesh::DofConstraints &constraints,
                                        const libMesh::DofMap &dof_map,
                                        const unsigned int variable_number,
                                        const libMesh::Elem * elem,
                                        const unsigned mesh_dim,
                                        std::vector<int> & index
                                        );

        static void compute_all_constraints(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map,
            libMesh::DofConstraints &constraints);

        static void compute_boundary_nodes(const libMesh::MeshBase &mesh, 
                                           libMesh::DofMap &dof_map,
                                           unsigned int sys_number, 
                                           unsigned int var_number,
                                           std::vector<int> & index);


        static void compute_boundary_nodes_top_parent(const libMesh::MeshBase &mesh, 
                                                       libMesh::DofMap &dof_map,
                                                       unsigned int sys_number, 
                                                       unsigned int var_number,
                                                       std::vector<int> & index);


        static void compute_boundary_nodes_to_skip(const libMesh::MeshBase &mesh, 
                                            libMesh::DofMap &dof_map,
                                            unsigned int sys_number, unsigned int var_number, 
                                            std::vector<int> & index);


        static void process_constraints (libMesh::MeshBase  &mesh, 
                                         libMesh::DofMap &dof_map, 
                                         libMesh::DofConstraints &_dof_constraints);

        static  void add_constraints_to_send_list(libMesh::DofMap &dof_map, 
                                                  libMesh::DofConstraints &_dof_constraints, std::vector<libMesh::dof_id_type> &_send_list );

        static void gather_constraints (libMesh::MeshBase  & mesh,
                                         std::set<libMesh::dof_id_type> & unexpanded_dofs, 
                                         libMesh::DofConstraints &_dof_constraints,
                                         libMesh::DofMap & dof_map,
                                         bool look_for_constrainees);

        static void allgather_recursive_constraints(libMesh::MeshBase  & mesh, 
                                                    libMesh::DofConstraints &_dof_constraints, 
                                                    libMesh::DofMap &dof_map);

        // static  void check_for_constraint_loops(const libMesh::DofConstraints & _dof_constraints, 
        //                                         libMesh::DofMap & dof_map);

        static void scatter_constraints(libMesh::MeshBase  & mesh, 
                                        libMesh::DofMap &dof_map, 
                                        libMesh::DofConstraints &_dof_constraints);



        static void merge_ghost_functor_outputs(libMesh::GhostingFunctor::map_type & elements_to_ghost,
                            std::set<libMesh::CouplingMatrix *> & temporary_coupling_matrices,
                            const std::set<libMesh::GhostingFunctor *>::iterator & gf_begin,
                            const std::set<libMesh::GhostingFunctor *>::iterator & gf_end,
                            const libMesh::MeshBase::const_element_iterator & elems_begin,
                            const libMesh::MeshBase::const_element_iterator & elems_end,
                            libMesh::processor_id_type p);
        static 
        void check_for_constraint_loops(libMesh::DofMap &dof_map, libMesh::DofConstraints &_dof_constraints);

        // static void compute_boundary_nodes(const libMesh::MeshBase &mesh, 
        //                                     const libMesh::DofMap &dof_map,
        //                                     unsigned int sys_number, unsigned int var_number, 
        //                                     std::vector<SizeType> & index)


        template<class ElementMatrix>
        static void construct_constraint_matrix(const libMesh::Elem *elem,
                                                const libMesh::DofMap &dof_map,
                                                const libMesh::DofConstraints &dof_constraints,
                                                ElementMatrix &mat,
                                                std::vector<libMesh::dof_id_type> &dof_indices,
                                                const bool called_recursively = false)
        {
            


            const std::size_t n_var = dof_map.n_variables();
            //dof_map.dof_indices(elem, dof_indices);

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
                    if (dof_map.is_constrained_dof(dof_indices[i])) 
                    {
                        auto pos = dof_constraints.find(dof_indices[i]);

                        if(pos == dof_constraints.end()) {
                            mat(i,i) = 1.;
                            continue;
                        }

                        const auto & constraint_row = pos->second;

                        for(const auto & item : constraint_row) {
                            const auto n_elem_dofs = static_cast<unsigned int>(dof_indices.size());

                           

                            for (unsigned int j=0; j != n_elem_dofs; j++) {
                                if (dof_indices[j] == item.first){
                                    mat(i, j) = 1.0 * item.second;
                                }
                            }
                        }
                    }

                    else 
                    {
                        mat(i,i) = 1.;
                    }
                }


                // ElementMatrix mat_new;

                // construct_constraint_matrix (elem, dof_map,  dof_constraints, mat_new, dof_indices, true);

                // if ((mat.n() == mat_new.m()) &&
                //       (mat_new.n() ==  dof_indices.size())) // If the constraint matrix
                //     mat.right_multiply(mat_new);           // is constrained...

                // libmesh_assert_equal_to (mat.n(),  dof_indices.size());
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

            //libMesh::DenseVector<double> tmp(vec);

            //tmp.zero();  

            mat.left_multiply_transpose(C);
            mat.right_multiply(C);

            //std::cout<<C<<std::endl;

            // mat.print_matlab();

            const uint ndofs = dof_indices.size();

            for(uint i = 0; i < ndofs; ++i) {
               
                if(dof_map.is_constrained_dof(dof_indices[i])) {
                    
                    auto it = dof_constraints.find(dof_indices[i]);
                    
                    if(it == dof_constraints.end()) continue;

                    for(uint j = 0; j < ndofs; j++) {
                        mat(i, j) = 0;
                    }

                    mat(i, i) = 1.0;

                    //tmp(i) = -1.0;

                    const auto &c_row = it->second;

                    for(const auto &item : c_row) {

                        for(uint j = 0; j < ndofs; j++) {
                            if(dof_indices[j] == item.first) {
                                //std::cout<<"dof_indices[i]:"<<dof_indices[i]<<" and hanging entry:"<<item.second<<std::endl;
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
