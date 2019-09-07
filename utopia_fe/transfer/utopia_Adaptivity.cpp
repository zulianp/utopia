/*  Created by Maria Nestola 07/09/2019*/

#include "utopia_Adaptivity.hpp"
#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"
#include "libmesh/fe_interface.h"

namespace utopia {

    void Adaptivity::compute_all_constraints(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        libMesh::DofConstraints &constraints)
    {
        using uint = unsigned int;
        auto el = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();
        
        for(; el != end_el; ++el) {
            const auto * elem = *el;
            auto * ele = *el;

            libmesh_assert (mesh.is_prepared());
            
            uint n_vars = dof_map.n_variables();
            for(uint var_num = 0; var_num < n_vars; ++var_num) {
                compute_constraints(constraints, dof_map, var_num, elem, mesh.mesh_dimension());
            }
        }

        std::cout << "--------------------------------------------------\n";
        std::cout<< "[Adaptivity::compute_all_constraints] n_constraints: " << constraints.size() << std::endl;
        std::cout << "--------------------------------------------------\n";
    }
    
    void Adaptivity::assemble_constraint(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num)
    {
    
        auto       el     = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();
        
        
        for ( ; el != end_el; ++el)
        {
            const auto * elem = *el;
            
            auto * ele = *el;
            libmesh_assert (mesh.is_prepared());
            

            for (unsigned int var_num=0; var_num<dof_map.n_variables();
                 ++var_num) {
                compute_constraints(dof_constraints_, dof_map,  var_num, elem, mesh.mesh_dimension());
            }
            
        }

        std::cout << "--------------------------------------------------\n";
        std::cout<< "[Adaptivity::assemble_constraint] n_constraints: " << dof_constraints_.size() << std::endl;
        std::cout << "--------------------------------------------------\n";
    }
    
    void Adaptivity::constraint_matrix(const LibMeshFunctionSpace &V, USparseMatrix &M, USparseMatrix &S)
    {
        const auto & mesh = V.mesh();
        unsigned int var_num = V.subspace_id();
        
        auto & dof_map = V.dof_map();
        constraint_matrix(mesh, dof_map, var_num, M, S);
    }

    void Adaptivity::constraint_matrix(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num, USparseMatrix &M, USparseMatrix &S)
    {

        assemble_constraint(mesh, dof_map, var_num);
        
        bool called_recursively = false;
    
        
        M = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);
        
        S = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);
        
        
        std::vector<libMesh::dof_id_type> elem_dofs;

        std::vector<libMesh::dof_id_type> I(1,0), J(1,0);

        std::vector<double> V(1, 0);
        
        auto      el     = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();
        
        
        libMesh::DenseMatrix<libMesh::Number> C;
        
        Write<USparseMatrix> w(M, utopia::GLOBAL_INSERT), w_2(S, utopia::GLOBAL_INSERT);
        
        for ( ; el != end_el; ++el)
        {
            const auto * elem = *el;
            
            auto * ele = *el;
            
            dof_map.dof_indices(elem, elem_dofs);
            
            std::set<libMesh::dof_id_type>  dof_set;
            
            bool we_have_constraints = false;

            for (const auto & dof : elem_dofs){
                
                if (dof_map.is_constrained_dof(dof))
                {
                    we_have_constraints = true;
                    
                    auto pos = dof_constraints_.find(dof);
                    
                    libmesh_assert (pos != dof_constraints_.end());
                    
                    const auto & constraint_row = pos->second;
                    
                    for (const auto & item : constraint_row)
                    dof_set.insert (item.first);
                    
                }
            }

            for (const auto & dof : elem_dofs)
            dof_set.erase (dof);


;
            
    
            if (!dof_set.empty() ||  
                !called_recursively) 
            {
                const unsigned int old_size =
                static_cast<unsigned int>(elem_dofs.size());
                
                elem_dofs.insert(elem_dofs.end(),
                                 dof_set.begin(), dof_set.end());
                
          
                C.resize (old_size,
                          static_cast<unsigned int>(elem_dofs.size()));
                
                for (unsigned int i=0; i != old_size; i++)
                if (dof_map.is_constrained_dof(elem_dofs[i]))
                {     
                    
                    I[0] = elem_dofs[i];
                    J[0] = elem_dofs[i];
                    V[0] = 0.0;

                    //S.set (elem_dofs[i],elem_dofs[i],0.0);
                    S.set_matrix(I, J, V);
                    
                    auto pos = dof_constraints_.find(elem_dofs[i]);
                    
                    libmesh_assert (pos != dof_constraints_.end());
                    
                    const auto & constraint_row = pos->second;
   
                    
                    for (const auto & item : constraint_row)
                    {
                        for (unsigned int j=0,
                             n_elem_dofs = static_cast<unsigned int>(elem_dofs.size());
                             j != n_elem_dofs; j++)
                        {
                            if (elem_dofs[j] == item.first){
                                C(i,j) = item.second;
                                  I[0] = elem_dofs[i];
                                  J[0] = elem_dofs[j];
                                  V[0] = item.second;
                                M.set_matrix(I,J,V);
                                //M.set (elem_dofs[i],elem_dofs[j],item.second);
                               
                                S.set_matrix(I,J,V);
                                //S.set (elem_dofs[i],elem_dofs[j],item.second);
                            }
                        }
                    }
                }
    
                else
                {
                    C(i,i) = 1.;
                    I[0] = elem_dofs[i];
                    J[0] = elem_dofs[i];
                    V[0] = 1.0;
                    M.set_matrix(I,J,V);
                    //M.set (elem_dofs[i],elem_dofs[i],1.0);
                }
            }
        }
    }
    

    void Adaptivity::compute_constraints(libMesh::DofConstraints &constraints,
                                         const libMesh::DofMap &dof_map,
                                         const unsigned int var_num,
                                         const libMesh::Elem * elem,
                                         const unsigned mesh_dim
                                         )
    {
        
        // Only constrain elements in 2,3D.
        if (mesh_dim == 1)
        return;
        
        // std::cout<<"lagrange_compute_constraints libmesh mio prima:"<<constraints.size()<<std::endl;
        libmesh_assert(elem);
        
        // Only constrain active and ancestor elements
        if (elem->subactive()) // if the element is subactive (i.e. has no active descendants)
        return;
        
        libMesh::FEType fe_type = dof_map.variable_type(var_num);
        fe_type.order = static_cast<libMesh::Order>(fe_type.order + elem->p_level());
        

        std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;
        std::unique_ptr<const libMesh::Elem> my_side, parent_side;

        for (auto s : elem->side_index_range())
        {
            if (elem->neighbor_ptr(s) != nullptr && elem->neighbor_ptr(s) != libMesh::remote_elem)
            {
                if (elem->neighbor_ptr(s)->level() < elem->level()) 
                
                {
                    const auto * parent = elem->parent();
                    libmesh_assert(parent);
                    
                    elem->build_side_ptr(my_side, s); 
                    
                    parent->build_side_ptr(parent_side, s);
                    
                    my_dof_indices.reserve (my_side->n_nodes()); 
                    
                    parent_dof_indices.reserve (parent_side->n_nodes());
                    
                    dof_map.dof_indices (my_side.get(), my_dof_indices,  var_num);

                    dof_map.dof_indices (parent_side.get(), parent_dof_indices, var_num);
                    
                    const unsigned int n_side_dofs = libMesh::FEInterface::n_dofs(mesh_dim-1, fe_type, my_side->type());
                    
                    const unsigned int n_parent_side_dofs = libMesh::FEInterface::n_dofs(mesh_dim-1, fe_type, parent_side->type());
                    

                    for (unsigned int my_dof=0; my_dof != n_side_dofs; my_dof++)
                    {
                        
                        libmesh_assert_less (my_dof, my_side->n_nodes());
                        
                        const libMesh::dof_id_type my_dof_g = my_dof_indices[my_dof];
                  
                        bool self_constraint = false;
                        for (unsigned int their_dof=0;
                             their_dof != n_parent_side_dofs; their_dof++)
                        {
                            libmesh_assert_less (their_dof, parent_side->n_nodes());
                            
                            const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];
                            
                            if (their_dof_g == my_dof_g)
                            {
                                self_constraint = true;
                                break;
                            }
                        }
                        
                        if (self_constraint) continue;
                        
                        libMesh::DofConstraintRow * constraint_row;
                     
                        constraint_row = &(constraints[my_dof_g]);

                        const libMesh::Point & support_point = my_side->point(my_dof);
                        
                        const libMesh::Point mapped_point = 
                        libMesh::FEInterface::inverse_map(mesh_dim-1, fe_type, parent_side.get(), support_point);
                        
                        for (unsigned int their_dof=0;
                             their_dof != n_parent_side_dofs; their_dof++)
                        {
                            libmesh_assert_less (their_dof, parent_side->n_nodes());
                            

                            const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];
                            
                            const libMesh::Real their_dof_value = libMesh::FEInterface::shape(mesh_dim-1,
                                                                            fe_type,
                                                                            parent_side->type(),
                                                                            their_dof,
                                                                            mapped_point);
                            if ((std::abs(their_dof_value) > 1.e-5) &&
                                (std::abs(their_dof_value) < .999))
                            {
                 
                              
                                
                                constraint_row->insert(std::make_pair (their_dof_g,
                                                                       their_dof_value));
                            }
                        }
                    }
                }
            }
        }
    }
}
     
