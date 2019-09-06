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
        // std::cout<<"Adaptivity::AssembleConstraint"<<std::endl;
        
        // Get a constant reference to the mesh object.
        
        //dof_map.createdof_constraints_(mesh, 0.0);
        
        
        auto       el     = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();
        
        // first we need to manually zero the matrix
        
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
        //LOG_SCOPE_IF("build_constraint_matrix()", "DofMap", !called_recursively);
        
        // std::cout<<"rintMatrix::constraint_matrix"<<dof_constraints_.size()<<std::endl;
        assemble_constraint(mesh, dof_map, var_num);
        
        // Create a set containing the DOFs we already depend on
        bool called_recursively = false;
        
       
        
        //const DofMap & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();
        
        // std::cout<<"begin"<<*dof_map.get_n_nz().begin()<<std::endl;
        
        // std::cout<<"end"<<*dof_map.get_n_nz().end()<<std::endl;
        
        // std::cout<<"size"<<dof_map.get_n_nz().size()<<std::endl;
        
        M = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);
        
        S = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);
        
        
        std::vector<libMesh::dof_id_type> elem_dofs;
        
        auto      el     = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();
        
        
        
        // ConstElemRange range (mesh.local_elements_begin(),
        //                      mesh.local_elements_end());
        
        // libMesh::Threads::parallel_for (range,
        //ComputeConstraints (dof_constraints_, dof_map,mesh,0);
        
        // first we need to manually zero the matrix
        
        
        libMesh::DenseMatrix<libMesh::Number> C;
        
        Write<USparseMatrix> w(M), w_2(S);
        
        for ( ; el != end_el; ++el)
        {
            const auto * elem = *el;
            
            auto * ele = *el;
            
            
            
            //ComputeConstraints(dof_constraints_, dof_map,  var_num, elem, 2);
            
            dof_map.dof_indices(elem, elem_dofs);
            
            typedef std::set<libMesh::dof_id_type> RCSet;
            
            RCSet dof_set;
            
            bool we_have_constraints = false;
            
            // std::cout<<"CREATE constraint_matrix"<<std::endl;
            
            // Next insert any other dofs the current dofs might be constrained
            // in terms of.  Note that in this case we may not be done: Those
            // may in turn depend on others.  So, we need to repeat this process
            // in that case until the system depends only on unconstrained
            // degrees of freedom.
            for (const auto & dof : elem_dofs){
                
                // std::cout<<"elem_dofs.size()"<<elem_dofs.size()<<"and constraint"<<dof_map.is_constrained_dof(dof)<<std::endl;
                
                if (dof_map.is_constrained_dof(dof))
                {
                    we_have_constraints = true;
                    
                    // If the DOF is constrained
                    
                    // std::cout<<"dof_map.is_constrained_dof(dof)"<<dof_map.is_constrained_dof(dof)<<std::endl;
                    
                    auto pos = dof_constraints_.find(dof);
                    
                    libmesh_assert (pos != dof_constraints_.end());
                    
                    const auto & constraint_row = pos->second;
                    
                    // Constraint rows in p refinement may be empty
                    //libmesh_assert (!constraint_row.empty());
                    
                    for (const auto & item : constraint_row)
                    dof_set.insert (item.first);
                    
                }
            }
            
            
            // std::cout<<"1::dof_map.is_constrained_dof(dof)"<<dof_map.is_constrained_dof(dof)<<std::endl;
            
            //   if (!we_have_constraints)
            //       return;
            
            
            // std::cout<<"2::dof_map.is_constrained_dof(dof)"<<dof_map.is_constrained_dof(dof)<<std::endl;
            
            for (const auto & dof : elem_dofs)
            dof_set.erase (dof);
            
            
            
            // If we added any DOFS then we need to do this recursively.
            // It is possible that we just added a DOF that is also
            // constrained!
            //
            // Also, we need to handle the special case of an element having DOFs
            // constrained in terms of other, local DOFs
            if (!dof_set.empty() ||  // case 1: constrained in terms of other DOFs
                !called_recursively) // case 2: constrained in terms of our own DOFs
            {
                const unsigned int old_size =
                static_cast<unsigned int>(elem_dofs.size());
                
                // Add new dependency dofs to the end of the current dof set
                elem_dofs.insert(elem_dofs.end(),
                                 dof_set.begin(), dof_set.end());
                
                // Now we can build the constraint matrix.
                // Note that resize also zeros for a libMesh::DenseMatrix<libMesh::Number>.
                C.resize (old_size,
                          static_cast<unsigned int>(elem_dofs.size()));
                
                // Create the C constraint matrix.
                for (unsigned int i=0; i != old_size; i++)
                if (dof_map.is_constrained_dof(elem_dofs[i]))
                {
                    // If the DOF is constrained
                    
                    S.set (elem_dofs[i],elem_dofs[i],0.0);
                    
                    auto pos = dof_constraints_.find(elem_dofs[i]);
                    
                    libmesh_assert (pos != dof_constraints_.end());
                    
                    const auto & constraint_row = pos->second;
                    
                    // p refinement creates empty constraint rows
                    //    libmesh_assert (!constraint_row.empty());
                    
                    for (const auto & item : constraint_row)
                    for (unsigned int j=0,
                         n_elem_dofs = static_cast<unsigned int>(elem_dofs.size());
                         j != n_elem_dofs; j++)
                    if (elem_dofs[j] == item.first){
                        C(i,j) = item.second;
                        M.set (elem_dofs[i],elem_dofs[j],item.second);
                        auto value = -1.0*item.second;
                        S.set (elem_dofs[i],elem_dofs[j],item.second);
                    }
                    
                }
                else
                {
                    C(i,i) = 1.;
                    M.set (elem_dofs[i],elem_dofs[i],1.0);
                }
            }
            // std::cout<<"current_elem_MOOSE: "<<ele[0]<<std::endl;
            
            // std::cout<<"C_MOOSE"<<"\n";
            // std::cout<<C<<std::endl;
            
            
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
        
        // Pull objects out of the loop to reduce heap operations
        std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;
        std::unique_ptr<const libMesh::Elem> my_side, parent_side;
        
        // Look at the element faces.  Check to see if we need to
        // build constraints.
        for (auto s : elem->side_index_range())
        if (elem->neighbor_ptr(s) != nullptr &&              //A const pointer to the $ i^{th} $ neighbor of this element
            elem->neighbor_ptr(s) != libMesh::remote_elem)
        if (elem->neighbor_ptr(s)->level() < elem->level()) // To check if the close element comes from refinement opeartion
        {
            // than this element.
            // Get pointers to the elements of interest and its parent.
            const auto * parent = elem->parent();
            
            //std::cout<<"parent->id()"<<parent->id()<<std::endl;
            
            // std::cout<<"elem->neighbor_ptr(s)->level()"<<elem->neighbor_ptr(s)->level()<<std::endl;
            
            // std::cout<<"elem->level()"<<elem->level()<<std::endl;
            
            // This can't happen...  Only level-0 elements have nullptr
            // parents, and no level-0 elements can be at a higher
            // level than their neighbors!
            libmesh_assert(parent);
            
            elem->build_side_ptr(my_side, s); // N.B.  calling build_side_ptr(s) on a 20-noded hex will build a 8-noded quadrilateral coincident with face s
            parent->build_side_ptr(parent_side, s);
            
            
            my_dof_indices.reserve (my_side->n_nodes()); //reallocation
            parent_dof_indices.reserve (parent_side->n_nodes());
            
            dof_map.dof_indices (my_side.get(), my_dof_indices,  // to get dofs of the faces on the element an the associated parent one
                                 var_num);
            dof_map.dof_indices (parent_side.get(), parent_dof_indices,
                                 var_num);
            
            const unsigned int n_side_dofs =
            libMesh::FEInterface::n_dofs(mesh_dim-1, fe_type, my_side->type());
            const unsigned int n_parent_side_dofs =
            libMesh::FEInterface::n_dofs(mesh_dim-1, fe_type, parent_side->type());
            for (unsigned int my_dof=0; my_dof != n_side_dofs; my_dof++)
            {
                libmesh_assert_less (my_dof, my_side->n_nodes());
                
                // My global dof index.
                const libMesh::dof_id_type my_dof_g = my_dof_indices[my_dof];
                
                // Hunt for "constraining against myself" cases before
                // we bother creating a constraint row
                bool self_constraint = false;
                for (unsigned int their_dof=0;
                     their_dof != n_parent_side_dofs; their_dof++)
                {
                    libmesh_assert_less (their_dof, parent_side->n_nodes());
                    
                    // Their global dof index.
                    const libMesh::dof_id_type their_dof_g =
                    parent_dof_indices[their_dof];
                    
                    if (their_dof_g == my_dof_g)
                    {
                        self_constraint = true;
                        break;
                    }
                }
                
                if (self_constraint)
                continue;
                
                libMesh::DofConstraintRow * constraint_row;
                
                // we may be running constraint methods concurrently
                // on multiple threads, so we need a lock to
                // ensure that this constraint is "ours"
                {
                    libMesh::Threads::spin_mutex::scoped_lock lock(libMesh::Threads::spin_mtx);
                    
                    // if (dof_map.is_constrained_dof(my_dof_g))
                    //  continue;
                    
                    //std::cout<<"dof_map.is_constrained_dof(my_dof_g)"<<dof_map.is_constrained_dof(my_dof_g)<<std::endl;
                    
                    constraint_row = &(constraints[my_dof_g]);
                    // libmesh_assert(constraint_row->empty());
                }
                
                // The support point of the DOF
                const libMesh::Point & support_point = my_side->point(my_dof);
                
                // Figure out where my node lies on their reference element.
                const libMesh::Point mapped_point = libMesh::FEInterface::inverse_map(mesh_dim-1, fe_type,
                                                                    parent_side.get(),
                                                                    support_point);
                
                // Compute the parent's side shape function values.
                for (unsigned int their_dof=0;
                     their_dof != n_parent_side_dofs; their_dof++)
                {
                    libmesh_assert_less (their_dof, parent_side->n_nodes());
                    
                    // Their global dof index.
                    const libMesh::dof_id_type their_dof_g =
                    parent_dof_indices[their_dof];
                    
                    const libMesh::Real their_dof_value = libMesh::FEInterface::shape(mesh_dim-1,
                                                                    fe_type,
                                                                    parent_side->type(),
                                                                    their_dof,
                                                                    mapped_point);
                    
                    // Only add non-zero and non-identity values
                    // for Lagrange basis functions.
                    if ((std::abs(their_dof_value) > 1.e-5) &&
                        (std::abs(their_dof_value) < .999))
                    {
                        
                        // std::cout<<"their_dof_g=>"<< their_dof_g <<std::endl;
                        // std::cout<<"their_dof_value=>"<< their_dof_value <<std::endl;
                        
                        // std::cout<<"elem->id()"<<elem->id()<<std::endl;
                        // std::cout<<"parent->id()"<<parent->id()<<std::endl;
                        
                        
                        constraint_row->insert(std::make_pair (their_dof_g,
                                                               their_dof_value));
                    }
#ifdef DEBUG
                    // Protect for the case u_i = 0.999 u_j,
                    // in which case i better equal j.
                    else if (their_dof_value >= .999)
                    libmesh_assert_equal_to (my_dof_g, their_dof_g);
#endif
                }
            }
        }
        
        // std::cout<<"lagrange_compute_constraints libmesh mio dopo:"<<constraints.size()<<std::endl;
        
        // dof_map.printdof_constraints_();
    } // lagrange_compute_constraints()
    
}
