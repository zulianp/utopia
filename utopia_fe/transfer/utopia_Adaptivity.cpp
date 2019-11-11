#include "utopia_Adaptivity.hpp"
#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_node.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/coupling_matrix.h"

namespace utopia {

    void Adaptivity::compute_all_constraints(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        libMesh::DofConstraints &constraints)
    {

      
        using uint = unsigned int;

        std::vector<int> index;

        libMesh::MeshBase &mesh_copy=const_cast<libMesh::MeshBase&>(mesh);

        libMesh::DofMap &dof_copy=const_cast<libMesh::DofMap&>(dof_map);

        uint vars = dof_map.n_variables();
        
        for(uint var_num = 0; var_num < vars; ++var_num) {

            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            fe_type.order = static_cast<libMesh::Order>(fe_type.order);

            if (fe_type.order>0){

                //compute_boundary_nodes(mesh, dof_copy, 0,0, index);
            }
        }
        
        auto el = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();
        
        for(; el != end_el; ++el) 
        {
            const auto * elem = *el;
            
            auto * ele = *el;

            libmesh_assert (mesh.is_prepared());
            
            uint n_vars = dof_map.n_variables();

            for(uint var_num = 0; var_num < n_vars; ++var_num) {

            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            fe_type.order = static_cast<libMesh::Order>(fe_type.order);

                if (fe_type.order>0){

                    compute_constraints(constraints, dof_map, var_num, elem, mesh.mesh_dimension(), index);
                }
            }
        }



        uint n_vars = dof_map.n_variables();

        for(uint var_num = 0; var_num < n_vars; ++var_num) 
        {

            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            fe_type.order = static_cast<libMesh::Order>(fe_type.order);

            if (fe_type.order>0)
            {
                cdof_copy.process_constraints(mesh_copy);
                process_constraints(mesh_copy, dof_copy, constraints);
            }
        }
      
        std::cout << "--------------------------------------------------\n";
        std::cout<< "[Adaptivity::compute_all_constraints] n_constraints: " << constraints.size() << std::endl;
        std::cout << "--------------------------------------------------\n";

      

        
    }
    
    void Adaptivity::assemble_constraint(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map)
    {
    
        auto       el     = mesh.active_local_elements_begin();
        
        const auto end_el = mesh.active_local_elements_end();

        dof_constraints_.clear();
        
        std::vector<int> index;

        libMesh::MeshBase &mesh_copy=const_cast<libMesh::MeshBase&>(mesh);

        libMesh::DofMap &dof_copy=const_cast<libMesh::DofMap&>(dof_map);

        //compute_boundary_nodes_to_skip(mesh, dof_copy, 0,0, index);

        libMesh::FEType fe_type = dof_map.variable_type(0);

        fe_type.order = static_cast<libMesh::Order>(fe_type.order);

        if(fe_type.order>0)
        {
        
        
            for ( ; el != end_el; ++el)
            {
                const auto * elem = *el;
                
                auto * ele = *el;
                libmesh_assert (mesh.is_prepared());
                

                for (unsigned int var_num=0; var_num<dof_map.n_variables();
                     ++var_num) {
                    compute_constraints(dof_constraints_, dof_map,  var_num, elem, mesh.mesh_dimension(), index);
                }
                
            }


        uint n_variables = dof_map.n_variables();

        for(uint var_num = 0; var_num < n_variables; ++var_num) 
        {

            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            fe_type.order = static_cast<libMesh::Order>(fe_type.order);

            if (fe_type.order>0)
            {
                cdof_copy.process_constraints(mesh_copy);

                Adaptivity::process_constraints(mesh_copy, dof_copy, dof_constraints_);
            }
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
       
        
        const auto & dof_map = V.dof_map();
        constraint_matrix(mesh, dof_map, M, S);

    }

    void Adaptivity::constraint_matrix(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, USparseMatrix &M, USparseMatrix &S)
    {
        using IndexArray = Traits<USparseMatrix>::IndexArray;
        
        std::vector<SizeType> index;

        assemble_constraint(mesh, dof_map);
        
        bool called_recursively = false;
    
        
        M = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);
        
        S = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);
        
        
        std::vector<libMesh::dof_id_type> elem_dofs;

        IndexArray I(1,0), J(1,0);

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

            for (const auto & dof : elem_dofs)
            {
                
                if (dof_map.is_constrained_dof(dof))
                {
                    we_have_constraints = true;
                    
                    auto pos = dof_constraints_.find(dof);
                    
                    if(pos == dof_constraints_.end()) continue;
                    
                    const auto & constraint_row = pos->second;
                    
                    for (const auto & item : constraint_row)
                    dof_set.insert (item.first);
                    
                }
            }

            for (const auto & dof : elem_dofs)
            dof_set.erase (dof);
            
    
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
                    
                    if(pos == dof_constraints_.end()) continue;
                    
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
                                         const unsigned mesh_dim,
                                         std::vector<int> & index
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
                        assert(my_dof < n_side_dofs);
                        
                        const libMesh::dof_id_type my_dof_g = my_dof_indices[my_dof];
                  
                        bool self_constraint = false;
                        
                        for (unsigned int their_dof=0; their_dof != n_parent_side_dofs; their_dof++)
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

                        //if (dof_map.is_constrained_dof(my_dof_g)) continue;

                        constraint_row = &(constraints[my_dof_g]);

                        const libMesh::Point & support_point = my_side->point(my_dof);
                        
                        const libMesh::Point mapped_point = 
                        libMesh::FEInterface::inverse_map(mesh_dim-1, fe_type, parent_side.get(), support_point);
                        
                        for (unsigned int their_dof=0;
                             their_dof != n_parent_side_dofs; their_dof++)
                        {
                            libmesh_assert_less (their_dof, parent_side->n_nodes());
                            

                            const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];
                            
                            const double their_dof_value = libMesh::FEInterface::shape(mesh_dim-1,
                                                                            fe_type,
                                                                            parent_side->type(),
                                                                            their_dof,
                                                                            mapped_point);
                            if ((std::abs(their_dof_value) > 1.e-5) &&
                                (std::abs(their_dof_value) < .999))
                            {
                                auto check = std::find(index.begin(), index.end(), their_dof_g) != index.end();
                 
                                if(!check)
                                constraint_row->insert(std::make_pair (their_dof_g,
                                                           their_dof_value));
                            
                      
                            }
                        }
                    }
                }
            }
        }
    }

    void Adaptivity::process_constraints (libMesh::MeshBase &mesh, libMesh::DofMap &dof_map, libMesh::DofConstraints &_dof_constraints)
    {

        using namespace libMesh;
        
        //dof_map.allgather_recursive_constraints(mesh);

        auto _primal_constraint_values = dof_map.get_primal_constraint_values();

        std::vector<int> index; 

        compute_boundary_nodes(mesh, dof_map, 0,0, index);

        typedef std::set<dof_id_type> RCSet;
        RCSet unexpanded_set;

        for (const auto & i : _dof_constraints)
        unexpanded_set.insert(i.first);

        while (!unexpanded_set.empty())
        for (RCSet::iterator i = unexpanded_set.begin();i != unexpanded_set.end(); /* nothing */)
        {
            // If the DOF is constrained
            DofConstraints::iterator pos = _dof_constraints.find(*i);

            libmesh_assert (pos != _dof_constraints.end());

            DofConstraintRow & constraint_row = pos->second;

            // DofConstraintValueMap::iterator rhsit =
            //   _primal_constraint_values.find(*i);
            // libMesh::Number constraint_rhs = (rhsit == _primal_constraint_values.end()) ?
            //   0 : rhsit->second;

            std::vector<dof_id_type> constraints_to_expand;

            for (const auto & item : constraint_row)
            {
                // if (item.first != *i && dof_map.is_constrained_dof(item.first))
                // {
                // unexpanded_set.insert(item.first);
                // constraints_to_expand.push_back(item.first);

                if (item.first != *i && dof_map.is_constrained_dof(item.first))
                {
                    bool check =true;

                    for(auto it=index.begin(); it < index.end(); ++it)
                    {
                        int b_id=*it;
                        
                        if(b_id==item.first) {

                            check = false;
                        }
                    }

                    if (check == true) {

                        unexpanded_set.insert(item.first);     

                        constraints_to_expand.push_back(item.first);
                    }

                }
            }

            for (const auto & expandable : constraints_to_expand)
            {
                const Real this_coef = constraint_row[expandable];

                DofConstraints::const_iterator
                subpos = _dof_constraints.find(expandable);

                libmesh_assert (subpos != _dof_constraints.end());

                if(subpos==_dof_constraints.end()) return;

                const DofConstraintRow & subconstraint_row = subpos->second;

                    for (const auto & item : subconstraint_row)
                    {
                        // Assert that the constraint does not form a cycle.
                        libmesh_assert(item.first != expandable);
                        constraint_row[item.first] += item.second * this_coef;
                    }

                // DofConstraintValueMap::const_iterator subrhsit =
                // _primal_constraint_values.find(expandable);
                // if (subrhsit != _primal_constraint_values.end())
                // constraint_rhs += subrhsit->second * this_coef;

                constraint_row.erase(expandable);
            }

            if (rhsit == _primal_constraint_values.end())
            {
                if (constraint_rhs != libMesh::Number(0))
                  _primal_constraint_values[*i] = constraint_rhs;
                else
                  _primal_constraint_values.erase(*i);
            }
            else
            {
                if (constraint_rhs != libMesh::Number(0))
                  rhsit->second = constraint_rhs;
                else
                  _primal_constraint_values.erase(rhsit);
            }

            if (constraints_to_expand.empty())
            i = unexpanded_set.erase(i);
            else
            ++i;
        }

        //dof_map.scatter_constraints(mesh);
        //dof_map.add_constraints_to_send_list();
    }
// <<<<<<< HEAD

//     void Adaptivity::allgather_recursive_constraints(libMesh::MeshBase & mesh, 
//                                                     libMesh::DofConstraints &_dof_constraints, 
//                                                     libMesh::DofMap &dof_map)
//     {
//        // std::cout<<"Adaptivity::allgather_recursive_constraints::BEGIN "<<std::endl;


//       if (dof_map.n_processors() == 1)
//         return;


//       unsigned int has_constraints = !_dof_constraints.empty();

//       dof_map.comm().max(has_constraints);

//       if (!has_constraints)
//         return;
//       {
//         std::map<libMesh::processor_id_type, std::set<libMesh::dof_id_type>> pushed_ids;



//         const unsigned int sys_num = dof_map.sys_number();

//         // Collect the constraints to push to each processor
//         for (auto & elem : as_range(mesh.active_not_local_elements_begin(),
//                                     mesh.active_not_local_elements_end()))
//           {
//             const unsigned short n_nodes = elem->n_nodes();

//             {
//               const unsigned int n_vars = elem->n_vars(sys_num);
//               for (unsigned int v=0; v != n_vars; ++v)
//                 {
//                   const unsigned int n_comp = elem->n_comp(sys_num,v);
//                   for (unsigned int c=0; c != n_comp; ++c)
//                     {
//                       const unsigned int id =
//                         elem->dof_number(sys_num,v,c);
//                       if (dof_map.is_constrained_dof(id))
//                         pushed_ids[elem->processor_id()].insert(id);
//                     }
//                 }
//             }

//             for (unsigned short n = 0; n != n_nodes; ++n)
//               {
//                 const libMesh::Node & node = elem->node_ref(n);
//                 const unsigned int n_vars = node.n_vars(sys_num);
//                 for (unsigned int v=0; v != n_vars; ++v)
//                   {
//                     const unsigned int n_comp = node.n_comp(sys_num,v);
//                     for (unsigned int c=0; c != n_comp; ++c)
//                       {
//                         const unsigned int id =
//                           node.dof_number(sys_num,v,c);
//                         if (dof_map.is_constrained_dof(id))
//                           pushed_ids[elem->processor_id()].insert(id);
//                       }
//                   }
//               }

//           }

//         // Rewrite those id sets as vectors for sending and receiving,
//         // then find the corresponding data for each id, then push it all.
//         std::map<libMesh::processor_id_type, std::vector<libMesh::dof_id_type>>
//           pushed_id_vecs, received_id_vecs;
//         for (auto & p : pushed_ids)
//           pushed_id_vecs[p.first].assign(p.second.begin(), p.second.end());

//         std::map<libMesh::processor_id_type, std::vector<std::vector<std::pair<libMesh::dof_id_type,double>>>>
//           pushed_keys_vals, received_keys_vals;

//         for (auto & p : pushed_id_vecs)
//         {
//             auto & keys_vals = pushed_keys_vals[p.first];
//             keys_vals.reserve(p.second.size());

//             for (auto & pushed_id : p.second)
//               {
//                 const libMesh::DofConstraintRow & row = _dof_constraints[pushed_id];
//                 keys_vals.emplace_back(row.begin(), row.end());

//               }
//         }

//         auto ids_action_functor =
//           [& received_id_vecs]
//           (libMesh::processor_id_type pid,
//            const std::vector<libMesh::dof_id_type> & data)
//           {
//             received_id_vecs[pid] = data;
//           };

//         libMesh::Parallel::push_parallel_vector_data
//           (dof_map.comm(), pushed_id_vecs, ids_action_functor);

//         auto keys_vals_action_functor =
//           [& received_keys_vals]
//           (libMesh::processor_id_type pid,
//            const std::vector<std::vector<std::pair<libMesh::dof_id_type,double>>> & data)
//           {
//             received_keys_vals[pid] = data;
//           };

//         libMesh::Parallel::push_parallel_vector_data
//           (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);



//         for (auto & p : received_id_vecs)
//           {
//             const libMesh::processor_id_type pid = p.first;
//             const auto & pushed_ids_to_me = p.second;
//             libmesh_assert(received_keys_vals.count(pid));

//             const auto & pushed_keys_vals_to_me = received_keys_vals.at(pid);


//             libmesh_assert_equal_to (pushed_ids_to_me.size(),
//                                      pushed_keys_vals_to_me.size());


//             for (auto i : libMesh::index_range(pushed_ids_to_me))
//               {
//                 libMesh::dof_id_type constrained = pushed_ids_to_me[i];

        
//                 if (!dof_map.is_constrained_dof(constrained))
//                   {
//                     libMesh::DofConstraintRow & row = _dof_constraints[constrained];
//                     for (auto & kv : pushed_keys_vals_to_me[i])
//                       {
//                         libmesh_assert_less(kv.first, dof_map.n_dofs());
//                         row[kv.first] = kv.second;
//                       }
//                   }
//               }
//           }
//       }


//       typedef std::set<libMesh::dof_id_type> DoF_RCSet;

//       DoF_RCSet unexpanded_dofs;

//       for (const auto & i : _dof_constraints) unexpanded_dofs.insert(i.first);


//       // std::cout<<"Adaptivity::allgather_recursive_constraints::END"<<std::endl;


//       gather_constraints(mesh, unexpanded_dofs, _dof_constraints, dof_map, false);

      
//     }



//     void Adaptivity::gather_constraints (libMesh::MeshBase & mesh,
//                                          std::set<libMesh::dof_id_type> & unexpanded_dofs, 
//                                          libMesh::DofConstraints & _dof_constraints,
//                                          libMesh::DofMap &dof_map,
//                                          bool look_for_constrainees)
//     {
      
//       //  std::cout<<"Adaptivity::gather_constraints::BEGIN "<<std::endl;  


//         typedef std::set<libMesh::dof_id_type> DoF_RCSet;

//         bool unexpanded_set_nonempty = !unexpanded_dofs.empty();

//         dof_map.comm().max(unexpanded_set_nonempty);

//         while (unexpanded_set_nonempty)
//         {
              
//               DoF_RCSet   dof_request_set;


//               std::map<libMesh::processor_id_type, std::vector<libMesh::dof_id_type>> requested_dof_ids;


//               std::map<libMesh::processor_id_type, libMesh::dof_id_type> dof_ids_on_proc;


//             for (const auto & unexpanded_dof : unexpanded_dofs)
//             {
//                   libMesh::DofConstraints::const_iterator
//                     pos = _dof_constraints.find(unexpanded_dof);

       
//                 if (pos == _dof_constraints.end())
//                 {
//                       if (!dof_map.local_index(unexpanded_dof) &&
//                           !_dof_constraints.count(unexpanded_dof) )
//                         dof_request_set.insert(unexpanded_dof);
//                 }

//                 else
//                 {
//                       const libMesh::DofConstraintRow & row = pos->second;
                     
//                     for (const auto & j : row)
//                     {
//                           const libMesh::dof_id_type constraining_dof = j.first;


//                           if (!dof_map.local_index(constraining_dof) &&
//                               !_dof_constraints.count(constraining_dof))
//                             dof_request_set.insert(constraining_dof);
//                     }
//                 }
//             }

//               unexpanded_dofs.clear();


//               libMesh::processor_id_type proc_id = 0;
//               for (const auto & i : dof_request_set)
//               {
//                   while (i >= dof_map.end_dof(proc_id))
//                   {
//                     proc_id++;
//                     dof_ids_on_proc[proc_id]++;
//                   }
//               }

//               for (auto & pair : dof_ids_on_proc)
//               {
//                   requested_dof_ids[pair.first].reserve(pair.second);
//               }

           
//               proc_id = 0;

//               for (const auto & i : dof_request_set)
//               {
//                   while (i >= dof_map.end_dof(proc_id))
//                   {
//                     proc_id++;
//                     requested_dof_ids[proc_id].push_back(i);
//                 }
//               }

//               unexpanded_set_nonempty = !unexpanded_dofs.empty();
//               dof_map.comm().max(unexpanded_set_nonempty);
//             }



//      //   std::cout<<"Adaptivity::gather_constraints::END "<<std::endl; 


//     }

//     void Adaptivity::add_constraints_to_send_list(libMesh::DofMap &dof_map, 
//                                                   libMesh::DofConstraints &_dof_constraints)
//     {
      
//         //std::cout<<"Adaptivity::add_constraints_to_send_list::BEGIN "<<std::endl; 


//         if (dof_map.n_processors() == 1) return;

          
//           unsigned int has_constraints = !_dof_constraints.empty();

//           dof_map.comm().max(has_constraints);

//         if (!has_constraints) return;

//         for (const auto & i : _dof_constraints)
//         {
//               libMesh::dof_id_type constrained_dof = i.first;

//               if (!dof_map.local_index(constrained_dof))
//                 continue;

//               const libMesh::DofConstraintRow & constraint_row = i.second;
              
//             for (const auto & j : constraint_row)
//             {
//                   libMesh::dof_id_type constraint_dependency = j.first;

//                   if (dof_map.local_index(constraint_dependency)) continue;

//                    std::vector<libMesh::dof_id_type> _send_list = dof_map.get_send_list();

//                   _send_list.push_back(constraint_dependency);
//             }
//         }       

//         //std::cout<<"Adaptivity::add_constraints_to_send_list::END "<<std::endl;  


//     }


//     void Adaptivity::scatter_constraints(libMesh::MeshBase & mesh, 
//                                          libMesh::DofMap &dof_map, 
//                                          libMesh::DofConstraints &_dof_constraints)
//     {




//       // std::cout<<"Adaptivity::scatter_constraints::BEGIN"<<std::endl;

//       // This function must be run on all processors at once
//       //parallel_object_only();

//       // Return immediately if there's nothing to gather
//       if (dof_map.n_processors() == 1)
//         return;

//       // We might get to return immediately if none of the processors
//       // found any constraints
//       unsigned int has_constraints = !_dof_constraints.empty();

        
//       dof_map.comm().max(has_constraints);

//       if (!has_constraints)
//         return;

//       // libMesh::Parallel::MessageTag range_tag = dof_map.comm().get_unique_tag();


//       std::map<libMesh::processor_id_type, std::set<libMesh::dof_id_type>> pushed_ids;

//       // Collect the dof constraints I need to push to each processor
//       libMesh::dof_id_type constrained_proc_id = 0;

//     for (auto & i : _dof_constraints)
//     {
//           const libMesh::dof_id_type constrained = i.first;
//           while (constrained >= dof_map.end_dof(constrained_proc_id))
//             constrained_proc_id++;

//           if (constrained_proc_id != dof_map.processor_id())
//             continue;

//           libMesh::DofConstraintRow & row = i.second;
//           for (auto & j : row)
//             {
//               const libMesh::dof_id_type constraining = j.first;

//               libMesh::processor_id_type constraining_proc_id = 0;
//               while (constraining >= dof_map.end_dof(constraining_proc_id))
//                 constraining_proc_id++;

//               if (constraining_proc_id != dof_map.processor_id() &&
//                   constraining_proc_id != constrained_proc_id)
//                 pushed_ids[constraining_proc_id].insert(constrained);
//             }
//         }

//       // Pack the dof constraint rows and rhs's to push

//       std::map<libMesh::processor_id_type,
//               std::vector<std::vector<std::pair<libMesh::dof_id_type, double>>>>
//         pushed_keys_vals, pushed_keys_vals_to_me;


//       auto keys_vals_action_functor =
//         [& pushed_keys_vals_to_me]
//         (libMesh::processor_id_type pid,
//          const std::vector<std::vector<std::pair<libMesh::dof_id_type, double>>> & data)
//         {
//           pushed_keys_vals_to_me[pid] = data;
//         };


//       libMesh::Parallel::push_parallel_vector_data
//         (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

      

//       typedef std::map<libMesh::dof_id_type, std::set<libMesh::dof_id_type>> DofConstrainsMap;
//       DofConstrainsMap dof_id_constrains;

//       for (auto & i : _dof_constraints)
//         {
//           const libMesh::dof_id_type constrained = i.first;
//           libMesh::DofConstraintRow & row = i.second;
//           for (const auto & j : row)
//             {
//               const libMesh::dof_id_type constraining = j.first;

//               libMesh::dof_id_type constraining_proc_id = 0;
//               while (constraining >= dof_map.end_dof(constraining_proc_id))
//                 constraining_proc_id++;

//               if (constraining_proc_id == dof_map.processor_id())
//                 dof_id_constrains[constraining].insert(constrained);
//             }
//         }

//       // Loop over all foreign elements, find any supporting our
//       // constrained dof indices.
//       pushed_ids.clear();

//       for (const auto & elem : as_range(mesh.active_not_local_elements_begin(),
//                                         mesh.active_not_local_elements_end()))
//         {
//           std::vector<libMesh::dof_id_type> my_dof_indices;
//           dof_map.dof_indices (elem, my_dof_indices);

//           for (const auto & dof : my_dof_indices)
//             {
//               DofConstrainsMap::const_iterator dcmi = dof_id_constrains.find(dof);
//               if (dcmi != dof_id_constrains.end())
//                 {
//                   for (const auto & constrained : dcmi->second)
//                     {
//                       libMesh::dof_id_type the_constrained_proc_id = 0;
//                       while (constrained >= dof_map.end_dof(the_constrained_proc_id))
//                         the_constrained_proc_id++;

//                       const libMesh::processor_id_type elemproc = elem->processor_id();
//                       if (elemproc != the_constrained_proc_id)
//                         pushed_ids[elemproc].insert(constrained);
//                     }
//                 }
//             }
//         }


//       pushed_keys_vals.clear();
//       pushed_keys_vals_to_me.clear();

//       libMesh::Parallel::push_parallel_vector_data
//         (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

//     }



// =======
    
// >>>>>>> remotes/origin/frac-refactor

    void Adaptivity::compute_boundary_nodes(const libMesh::MeshBase &mesh, 
                                            libMesh::DofMap &dof_map,
                                            unsigned int sys_number, unsigned int var_number, 
                                            std::vector<int> & index)
    {

       
       // std::cout<<"Adaptivity::compute_boundary_nodes::Begin "<<std::endl; 
       
       auto on_boundary = libMesh::MeshTools::find_boundary_nodes(mesh);     

       std::vector<int> dirichlet_id, index_local;

       index_local.clear();

       dirichlet_id.clear();

       // dirichlet_id.push_back(2);

       // dirichlet_id.push_back(4);

       index.clear(); 

       if(mesh.mesh_dimension()<3)
       {

            libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();
            
            const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();
            
            for ( ; it != end_it; ++it)
            {
                const libMesh::Elem * ele = *it;

                for(int kk=0; kk<ele->n_sides(); kk++) {       
                    
                    auto neigh = ele->neighbor_ptr(kk);    

                    if (neigh != libMesh::remote_elem && mesh.get_boundary_info().boundary_id(ele, kk)>0)
                    {
                        auto side = ele->build_side_ptr(kk);

                        index_local.clear();

                        for (int ll=0; ll<ele->n_nodes(); ll++)
                        {

                           const libMesh::Node * node = ele->node_ptr(ll);

                           const libMesh::dof_id_type node_dof = node->dof_number(sys_number, var_number, 0);                

                            if(on_boundary.count(node->id()) && dof_map.is_constrained_dof(node_dof)) 
                            {
                                   
                                index_local.push_back(node_dof);
           
                            }

                        }

                        if(index_local.size()==side->n_nodes())
                        {

                           index.insert(index.end(), index_local.begin(), index_local.end());
                        }
                    }
                }
            }
        }

       else
       {
            {
                libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

                libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();
              
                const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();

                std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;
            
                std::unique_ptr<const libMesh::Elem> my_side, parent_side_0;

                libMesh::FEType fe_type = dof_map.variable_type(0);
          

              
                for ( ; it != end_it; ++it)
                {
                    const libMesh::Elem * ele = *it; 

                    const auto *ele_parent_0 = ele->top_parent();

                    for(int jj=0; jj<ele_parent_0->n_sides(); jj++) 
                    {

                                
                      libmesh_assert(ele_parent_0);

                      auto parent_side_0 = ele_parent_0->build_side_ptr(jj);

                      index_local.clear();

                       for (int ll=0; ll<parent_side_0->n_nodes(); ll++)
                       {
                    
                            const libMesh::Node * node_0 = parent_side_0->node_ptr(ll);

                            const libMesh::dof_id_type node_dof_0 = node_0->dof_number(sys_number, var_number, 0); 
                           
                            if(dof_map.is_constrained_dof(node_dof_0)) {
                                
                                index_local.push_back(node_dof_0);

                                auto valpos = rhs_values.find(node_dof_0);

                                index.push_back(node_dof_0);                                               
                            }
                        }

                           if(index_local.size()==parent_side_0->n_nodes()){

                            auto bc_id = mesh.get_boundary_info().boundary_id(ele_parent_0,jj);

                            auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());


                            if (!check) dirichlet_id.push_back(bc_id);
                        }
                    }
                }
            }

            libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();
                
            const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();
                
            for ( ; it != end_it; ++it)
            {
                const libMesh::Elem * ele = *it;

                for(int kk=0; kk<ele->n_sides(); kk++) 
                {     

                    auto neigh = ele->neighbor_ptr(kk); 

                    auto bc_id = mesh.get_boundary_info().boundary_id(ele,kk);

                    auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());

                    if(check)
                    {
                         index_local.clear();

                         auto side = ele->build_side_ptr(kk);

                        for (int ll=0; ll<side->n_nodes(); ll++)
                        {
                          
                            const libMesh::Node * node = side->node_ptr(ll);

                            const libMesh::dof_id_type node_dof = node->dof_number(sys_number, var_number, 0); 

                            if(on_boundary.count(node->id()) && dof_map.is_constrained_dof(node_dof)) index.push_back(node_dof);
                        }
                    }
                }
            }
        }


    // std::cout<<"Adaptivity::compute_boundary_nodes::END "<<std::endl; 
    }


    void Adaptivity::compute_boundary_nodes_to_skip(const libMesh::MeshBase &mesh, 
                                            libMesh::DofMap &dof_map,
                                            unsigned int sys_number, unsigned int var_number, 
                                            std::vector<int> & index)
    {

        // std::cout<<"Adaptivity::compute_boundary_nodes_to_skip::BEGIN "<<std::endl;

       // Only constrain elements in 2,3D.
        if (mesh.mesh_dimension() == 1)
        return;

        unsigned int mesh_dim = mesh.mesh_dimension();
        
        // std::cout<<"lagrange_compute_constraints libmesh mio prima:"<<constraints.size()<<std::endl;
        //libmesh_assert(elem);

        std::vector<int> index_self;

        index_self.clear();

        unsigned int var_num = 0;


        libMesh::MeshBase::const_element_iterator it_0 = mesh.active_elements_begin();
      
        const libMesh::MeshBase::const_element_iterator end_it_0 = mesh.active_elements_end();        
        


        for ( ; it_0 != end_it_0; ++it_0)
        {
            const libMesh::Elem * elem = *it_0;

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
                        assert(my_dof < n_side_dofs);
                        
                        const libMesh::dof_id_type my_dof_g = my_dof_indices[my_dof];
                  
                        bool self_constraint = false;
                        
                            for (unsigned int their_dof=0; their_dof != n_parent_side_dofs; their_dof++)
                            {
                                libmesh_assert_less (their_dof, parent_side->n_nodes());
                                
                                const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];
                                
                                if (their_dof_g == my_dof_g)
                                {
                                    //self_constraint = true;

                                    auto check_2 = (std::find(index_self.begin(), index_self.end(), their_dof_g) != index_self.end());

                                    if(!check_2) {index_self.push_back(their_dof_g);

                                    //std::cout<<"their_dof_g"<<their_dof_g<<std::endl;

                                    }

                                    //index_self.push_back(their_dof_g);

                                    //std::cout<<"their_dof_g"<<their_dof_g<<std::endl;
                                    //break;
                                }
                            }                        
                        }
                    }
                }
            }
        }
            // std::cout<<"Adaptivity::compute_boundary_nodes::Begin "<<std::endl; 
       
        auto on_boundary = libMesh::MeshTools::find_boundary_nodes(mesh);     

        std::vector<int> dirichlet_id, index_local, tmp;

        index_local.clear();

        dirichlet_id.clear();

        index.clear(); 

            


        libMesh::MeshBase::const_element_iterator it_1 = mesh.active_elements_begin();
      
        const libMesh::MeshBase::const_element_iterator end_it_1 = mesh.active_elements_end();
    
        std::unique_ptr<const libMesh::Elem> parent_side_0_new;
          
        for ( ; it_1 != end_it_1; ++it_1)
        {
            const libMesh::Elem * ele = *it_1; 

            const auto *ele_parent_0 = ele->top_parent();

            for(int jj=0; jj<ele_parent_0->n_sides(); jj++) 
            {

                        
              libmesh_assert(ele_parent_0);

              auto parent_side_0_new = ele_parent_0->build_side_ptr(jj);

              index_local.clear();

               for (int ll=0; ll<parent_side_0_new->n_nodes(); ll++)
               {
            
                    const libMesh::Node * node_0 = parent_side_0_new->node_ptr(ll);

                    const libMesh::dof_id_type node_dof_0 = node_0->dof_number(sys_number, var_number, 0); 
                   
                    if(dof_map.is_constrained_dof(node_dof_0)) {
                        
                        index_local.push_back(node_dof_0);                                            
                    }
                }

                   if(index_local.size()==parent_side_0_new->n_nodes()){

                    auto bc_id = mesh.get_boundary_info().boundary_id(ele_parent_0,jj);

                    auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());


                    if (!check) dirichlet_id.push_back(bc_id);
                }
            }
        }
        

        libMesh::MeshBase::const_element_iterator it_2 = mesh.active_elements_begin();
            
        const libMesh::MeshBase::const_element_iterator end_it_2 = mesh.active_elements_end();
            
        for ( ; it_2 != end_it_2; ++it_2)
        {
            const libMesh::Elem * ele = *it_2;

           // if (ele->subactive()) 
            {// if the element is subactive (i.e. has no active descendants)
 
                for(int kk=0; kk<ele->n_sides(); kk++) 
                {     

                    auto neigh = ele->neighbor_ptr(kk);
                   if(ele->neighbor_ptr(kk) == nullptr && ele->neighbor_ptr(kk) !=  libMesh::remote_elem)
                   {

                        auto bc_id = mesh.get_boundary_info().boundary_id(ele,kk);

                        auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());

                        if(check)
                        {
                             index_local.clear();

                             auto side = ele->build_side_ptr(kk); 

                            for (int ll=0; ll<side->n_nodes(); ll++)
                            {
                              
                                const libMesh::Node * node = side->node_ptr(ll);

                                const libMesh::dof_id_type node_dof = node->dof_number(sys_number, var_number, 0); 

                                if(on_boundary.count(node->id()) && dof_map.is_constrained_dof(node_dof)) {

                                  auto check_2 = (std::find(tmp.begin(), tmp.end(), node_dof) != tmp.end());

                                   if(!check_2) {tmp.push_back(node_dof);}
                                    //std::cout<<"tmp"<<node_dof<<std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }


        for(auto it=tmp.begin(); it!=tmp.end(); ++it){

            auto check = (std::find(index_self.begin(), index_self.end(), *it) != index_self.end());

            if(!check){

                auto check_2 = (std::find(index.begin(), index.end(), *it) != index.end());

                if(!check_2) {
                    index.push_back(*it);
                    //std::cout<<"index_to_skip"<<*it<<std::endl;
                }
            }
        }
    }



    // void Adaptivity::process_constraints (libMesh::MeshBase &mesh, libMesh::DofMap &dof_map, libMesh::DofConstraints &_dof_constraints)
    // {

    //     std::cout<<"Adaptivity::process_constraints::BEGIN "<<std::endl;

       
    //     std::vector<int> index; 

    //     std::vector<int> index_2; 

    //     using namespace libMesh;

    //     //check_for_constraint_loops(dof_map, _dof_constraints);

    //     dof_map.allgather_recursive_constraints(mesh);

    //     //compute_boundary_nodes_to_skip(mesh, dof_map, 0,0, index);

    //     compute_boundary_nodes(mesh, dof_map, 0,0, index);         
    
    //     typedef std::set<dof_id_type> RCSet;

    //     std::set<dof_id_type> unexpanded_set;

    //     for (const auto & i : _dof_constraints)
    //          unexpanded_set.insert(i.first);

    //     while (!unexpanded_set.empty())
    //         for (RCSet::iterator i = unexpanded_set.begin(); i != unexpanded_set.end();)
    //         {
           
    //             DofConstraints::iterator pos = _dof_constraints.find(*i);

    //             DofConstraintRow & constraint_row = pos->second;

    //             std::vector<dof_id_type> constraints_to_expand;

    //             for (const auto & item : constraint_row)
    //             {
    //                 if (item.first != *i && dof_map.is_constrained_dof(item.first))
    //                 {
    //                     bool check =true;

    //                     for(auto it=index.begin(); it < index.end(); ++it)
    //                     {
    //                         int b_id=*it;
                            
    //                         if(b_id==item.first) {

    //                             check = false;
    //                         }
    //                     }
             
    //                     if (check == true) {

    //                         unexpanded_set.insert(item.first);     

    //                         constraints_to_expand.push_back(item.first);
    //                     }
                    
    //                 }
    //             }

                

    //             for (const auto & expandable : constraints_to_expand)
    //             {
    //                 const Real this_coef = constraint_row[expandable];

    //                 DofConstraints::const_iterator
    //                                 subpos = _dof_constraints.find(expandable);


    //                 const DofConstraintRow & subconstraint_row = subpos->second;
                   
    //                 libMesh::processor_id_type constraining_proc_id = 0;
                    
    //                 for (const auto & item : subconstraint_row)
    //                 {
    //                    // Assert that the constraint does not form a cycle.
    //                    //if(item.first == expandable) return; 

                        
                        
    //                     // const dof_id_type constraining = item.first;
                       
    //                     // while (constraining >= dof_map.end_dof(constraining_proc_id))
    //                     //   constraining_proc_id++;

    //                     // if (constraining_proc_id != dof_map.processor_id())
    //                     // continue;

    //                     constraint_row[item.first] += item.second * this_coef;
    //                 }

    //                 // DofConstraintValueMap::const_iterator subrhsit =
    //                 //   _primal_constraint_values.find(expandable);
    //                 // if (subrhsit != _primal_constraint_values.end())
    //                 //   constraint_rhs += subrhsit->second * & dof_map_coef;

    //                 constraint_row.erase(expandable);
    //             }

    //             // if (rhsit == _primal_constraint_values.end())
    //             //   {
    //             //     if (constraint_rhs != Number(0))
    //             //       _primal_constraint_values[*i] = constraint_rhs;
    //             //     else
    //             //       _primal_constraint_values.erase(*i);
    //             //   }
    //             // else
    //             //   {
    //             //     if (constraint_rhs != Number(0))
    //             //       rhsit->second = constraint_rhs;
    //             //     else
    //             //       _primal_constraint_values.erase(rhsit);
    //             //   }

    //             if (constraints_to_expand.empty()) i = unexpanded_set.erase(i);
    //             else ++i;
    //         }
        


    //     //dof_map.prepare_send_list();


    //     // dof_map.reinit_send_list(mesh);

    //     // std::shared_ptr<std::vector<libMesh::dof_id_type>> _send_list = std::make_shared<std::vector<libMesh::dof_id_type>>(dof_map.get_send_list());

    //     dof_map.scatter_constraints(mesh);
      
    //     dof_map.add_constraints_to_send_list();
 
        
    //     std::cout<<"Adaptivity::process_constraints::END "<<std::endl;    

    // }
 

    
   //  void Adaptivity::allgather_recursive_constraints(libMesh::MeshBase & mesh, 
   //                                                  libMesh::DofConstraints &_dof_constraints, 
   //                                                  libMesh::DofMap &dof_map)
   //  {
       

   //      // std::cout<<"Adaptivity::allgather_recursive_constraints::BEGIN "<<std::endl;
   //      // parallel_object_only();

   //      using namespace libMesh;

   //      auto _primal_constraint_values = dof_map.get_primal_constraint_values();

   //      //& dof_map function must be run on all processors at once
   //      //parallel_object_only();

   //      std::cout<<"I am in DofMap::allgather_recursive_constraints"<<std::endl;

   //      // Return immediately if there's nothing to gather
   //      if (dof_map.n_processors() == 1)
   //      return;

   //      // We might get to return immediately if none of the processors
   //      // found any constraints
   //      unsigned int has_constraints = !_dof_constraints.empty();
        
   //      dof_map.comm().max(has_constraints);
        
   //      if (!has_constraints)
   //      return;

   //      // If we have heterogenous adjoint constraints we need to
   //      // communicate those too.
   //      const unsigned int max_qoi_num = 0;


   //      // We might have calculated constraints for constrained dofs
   //      // which have support on other processors.
   //      // Push these out first.
   //      {
   //      std::map<processor_id_type, std::set<dof_id_type>> pushed_ids;

   //      const unsigned int sys_num = dof_map.sys_number();

   //      // Collect the constraints to push to each processor
   //      for (auto & elem : as_range(mesh.active_not_local_elements_begin(),
   //                                  mesh.active_not_local_elements_end()))
   //        {
   //          const unsigned short n_nodes = elem->n_nodes();

   //          // Just checking dof_indices on the foreign element isn't
   //          // enough.  Consider a central hanging node between a coarse
   //          // Q2/Q1 element and its finer neighbors on a higher-ranked
   //          // processor.  The coarse element's processor will own the node,
   //          // and will thereby own the pressure dof on that node, despite
   //          // the fact that that pressure dof doesn't directly exist on the
   //          // coarse element!
   //          //
   //          // So, we loop through dofs manually.

   //          {
   //            const unsigned int n_vars = elem->n_vars(sys_num);
   //            for (unsigned int v=0; v != n_vars; ++v)
   //              {
   //                const unsigned int n_comp = elem->n_comp(sys_num,v);
   //                for (unsigned int c=0; c != n_comp; ++c)
   //                  {
   //                    const unsigned int id =
   //                      elem->dof_number(sys_num,v,c);
   //                    if (dof_map.is_constrained_dof(id))
   //                      pushed_ids[elem->processor_id()].insert(id);
   //                  }
   //              }
   //          }

   //          for (unsigned short n = 0; n != n_nodes; ++n)
   //            {
   //              const Node & node = elem->node_ref(n);
   //              const unsigned int n_vars = node.n_vars(sys_num);
   //              for (unsigned int v=0; v != n_vars; ++v)
   //                {
   //                  const unsigned int n_comp = node.n_comp(sys_num,v);
   //                  for (unsigned int c=0; c != n_comp; ++c)
   //                    {
   //                      const unsigned int id =
   //                        node.dof_number(sys_num,v,c);
   //                      if (dof_map.is_constrained_dof(id))
   //                        pushed_ids[elem->processor_id()].insert(id);
   //                    }
   //                }
   //            }
   //        }

   //      // Rewrite those id sets as vectors for sending and receiving,
   //      // then find the corresponding data for each id, then push it all.
   //      std::map<processor_id_type, std::vector<dof_id_type>> pushed_id_vecs, received_id_vecs;

   //      for (auto & p : pushed_ids) pushed_id_vecs[p.first].assign(p.second.begin(), p.second.end());

   //      std::map<processor_id_type, std::vector<std::vector<std::pair<dof_id_type,Real>>>> pushed_keys_vals, received_keys_vals;

   //      std::map<processor_id_type, std::vector<std::vector<libMesh::Number>>> pushed_rhss, received_rhss;

   //      for (auto & p : pushed_id_vecs)
   //      {
   //          auto & keys_vals = pushed_keys_vals[p.first];

   //          keys_vals.reserve(p.second.size());

   //          auto & rhss = pushed_rhss[p.first];

   //          rhss.reserve(p.second.size());

   //          for (auto & pushed_id : p.second)
   //          {
   //              const DofConstraintRow & row = _dof_constraints[pushed_id];

   //              keys_vals.emplace_back(row.begin(), row.end());

   //              rhss.push_back(std::vector<libMesh::Number>(max_qoi_num+1));

   //              std::vector<libMesh::Number> & rhs = rhss.back();

   //              DofConstraintValueMap::const_iterator rhsit =
   //                _primal_constraint_values.find(pushed_id);
   //              rhs[max_qoi_num] =
   //                (rhsit == _primal_constraint_values.end()) ?
   //                0 : rhsit->second;
   //              // for (unsigned int q = 0; q != max_qoi_num; ++q)
   //              //   {
   //              //     AdjointDofConstraintValues::const_iterator adjoint_map_it =
   //              //       _adjoint_constraint_values.find(q);

   //              //     if (adjoint_map_it == _adjoint_constraint_values.end())
   //              //       continue;

   //              //     const DofConstraintValueMap & constraint_map =
   //              //       adjoint_map_it->second;

   //              //     DofConstraintValueMap::const_iterator adj_rhsit =
   //              //       constraint_map.find(pushed_id);

   //              //     rhs[q] =
   //              //       (adj_rhsit == constraint_map.end()) ?
   //              //       0 : adj_rhsit->second;
   //              //   }
   //            }
   //        }

   //      auto ids_action_functor =
   //        [& received_id_vecs]
   //        (processor_id_type pid,
   //         const std::vector<dof_id_type> & data)
   //        {
   //          received_id_vecs[pid] = data;
   //        };

   //      Parallel::push_parallel_vector_data(dof_map.comm(), pushed_id_vecs, ids_action_functor);

   //      auto keys_vals_action_functor =
   //        [& received_keys_vals]
   //        (processor_id_type pid,
   //         const std::vector<std::vector<std::pair<dof_id_type,Real>>> & data)
   //        {
   //          received_keys_vals[pid] = data;
   //        };

   //      Parallel::push_parallel_vector_data
   //        (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

   //      auto rhss_action_functor =
   //        [& received_rhss]
   //        (processor_id_type pid,
   //         const std::vector<std::vector<libMesh::Number>> & data)
   //        {
   //          received_rhss[pid] = data;
   //        };

   //      Parallel::push_parallel_vector_data(dof_map.comm(), pushed_rhss, rhss_action_functor);

   //      // Now we have all the DofConstraint rows and rhs values received
   //      // from others, so add the DoF constraints that we've been sent


   //      // Add all the dof constraints that I've been sent
   //      for (auto & p : received_id_vecs)
   //        {
   //          const processor_id_type pid = p.first;

   //          const auto & pushed_ids_to_me = p.second;

   //          libmesh_assert(received_keys_vals.count(pid));

   //          libmesh_assert(received_rhss.count(pid));

   //          const auto & pushed_keys_vals_to_me = received_keys_vals.at(pid);

   //          const auto & pushed_rhss_to_me = received_rhss.at(pid);

   //          libmesh_assert_equal_to (pushed_ids_to_me.size(),
   //                                   pushed_keys_vals_to_me.size());
   //          libmesh_assert_equal_to (pushed_ids_to_me.size(),
   //                                   pushed_rhss_to_me.size());

   //          for (auto i : index_range(pushed_ids_to_me))
   //          {
   //              dof_id_type constrained = pushed_ids_to_me[i];

   //              // If we don't already have a constraint for & dof_map dof,
   //              // add the one we were sent
   //              if (!dof_map.is_constrained_dof(constrained))
   //              {
   //                  DofConstraintRow & row = _dof_constraints[constrained];

   //                  for (auto & kv : pushed_keys_vals_to_me[i])
   //                  {
   //                      libmesh_assert_less(kv.first, dof_map.n_dofs());
   //                      row[kv.first] = kv.second;
   //                  }

   //                  const libMesh::Number primal_rhs = pushed_rhss_to_me[i][max_qoi_num];

   //                  if (libmesh_isnan(primal_rhs)) libmesh_assert(pushed_keys_vals_to_me[i].empty());

   //                  if (primal_rhs != libMesh::Number(0))
   //                    _primal_constraint_values[constrained] = primal_rhs;
   //                  else
   //                    _primal_constraint_values.erase(constrained);

   //                  // for (unsigned int q = 0; q != max_qoi_num; ++q)
   //                  // {
   //                      //     AdjointDofConstraintValues::iterator adjoint_map_it =
   //                      //       _adjoint_constraint_values.find(q);

   //                      //     const Number adj_rhs = pushed_rhss_to_me[i][q];

   //                      //     if ((adjoint_map_it == _adjoint_constraint_values.end()) &&
   //                      //         adj_rhs == Number(0))
   //                      //       continue;

   //                      //     if (adjoint_map_it == _adjoint_constraint_values.end())
   //                      //       adjoint_map_it = _adjoint_constraint_values.insert
   //                      //         (std::make_pair(q,DofConstraintValueMap())).first;

   //                      //     DofConstraintValueMap & constraint_map =
   //                      //       adjoint_map_it->second;

   //                      //     if (adj_rhs != Number(0))
   //                      //       constraint_map[constrained] = adj_rhs;
   //                      //     else
   //                      //       constraint_map.erase(constrained);
   //                  //   }
   //                  }
   //              }
   //          }
   //      }

   //      // Now start checking for any other constraints we need
   //      // to know about, requesting them recursively.

   //      // Create sets containing the DOFs and nodes we already depend on
   //      typedef std::set<dof_id_type> DoF_RCSet;
       
   //      DoF_RCSet unexpanded_dofs;

   //      for (const auto & i : _dof_constraints) unexpanded_dofs.insert(i.first);

   //      // Gather all the dof constraints we need
   //      gather_constraints(mesh, unexpanded_dofs, _dof_constraints, dof_map, false);

       
   //  }

   //  void Adaptivity::gather_constraints (libMesh::MeshBase  & mesh,
   //                                   std::set<libMesh::dof_id_type> & unexpanded_dofs, 
   //                                   libMesh::DofConstraints &_dof_constraints,
   //                                   libMesh::DofMap & dof_map,
   //                                   bool look_for_constrainees)
   //  {
        
   //      using namespace libMesh;

   //      auto _primal_constraint_values = dof_map.get_primal_constraint_values();

   //      typedef std::set<dof_id_type> DoF_RCSet;

   //      const unsigned int max_qoi_num = 0;
    
   //      bool unexpanded_set_nonempty = !unexpanded_dofs.empty();
    
   //      dof_map.comm().max(unexpanded_set_nonempty);

   //      while (unexpanded_set_nonempty)
   //      {
            
   //          DoF_RCSet   dof_request_set;

   //          std::map<processor_id_type, std::vector<dof_id_type>>
   //          requested_dof_ids;

   //          std::map<processor_id_type, dof_id_type>
   //          dof_ids_on_proc;

   //          for (const auto & unexpanded_dof : unexpanded_dofs)
   //          {
   //              DofConstraints::const_iterator pos = _dof_constraints.find(unexpanded_dof);

            
   //              if (pos == _dof_constraints.end())
   //              {
   //                  if (!dof_map.local_index(unexpanded_dof) &&
   //                      !_dof_constraints.count(unexpanded_dof) )
   //                  dof_request_set.insert(unexpanded_dof);
   //              }
   //              else
   //              {
   //                  const DofConstraintRow & row = pos->second;
   //                  for (const auto & j : row)
   //                  {
   //                    const dof_id_type constraining_dof = j.first;

   //                    // If it's non-local and we haven't already got a
   //                    // constraint for it, we might need to ask for one
   //                    if (!dof_map.local_index(constraining_dof) &&
   //                        !_dof_constraints.count(constraining_dof))
   //                      dof_request_set.insert(constraining_dof);
   //                  }
   //              }
   //          }

           
   //          unexpanded_dofs.clear();

   //          processor_id_type proc_id = 0;
   //          for (const auto & i : dof_request_set)
   //          {
   //              while (i >= dof_map.end_dof(proc_id)) proc_id++;
   //              dof_ids_on_proc[proc_id]++;
   //          }

   //          for (auto & pair : dof_ids_on_proc)
   //          {
   //            requested_dof_ids[pair.first].reserve(pair.second);
   //          }

   //          // Prepare each processor's request set
   //          proc_id = 0;

   //          for (const auto & i : dof_request_set)
   //          {
   //             while (i >= dof_map.end_dof(proc_id)) proc_id++;
               
   //             requested_dof_ids[proc_id].push_back(i);
   //          }

   //          typedef std::vector<std::pair<dof_id_type, Real>> row_datum;

   //          typedef std::vector<libMesh::Number> rhss_datum;

   //          auto row_gather_functor =
   //          [&dof_map,
   //          &_dof_constraints]
   //          (processor_id_type,
   //           const std::vector<dof_id_type> & ids,
   //           std::vector<row_datum> & data)
   //          {
   //              const std::size_t query_size = ids.size();

   //              data.resize(query_size);
             
   //              for (std::size_t i=0; i != query_size; ++i)
   //              {
   //                  dof_id_type constrained = ids[i];

   //                  if (_dof_constraints.count(constrained))
   //                  {
   //                      DofConstraintRow & row = _dof_constraints[constrained];
                        
   //                      std::size_t row_size = row.size();
                        
   //                      data[i].reserve(row_size);
                        
   //                      for (const auto & j : row)
   //                      {
   //                        data[i].push_back(j);

   //                        libmesh_assert(j.first != DofObject::invalid_id);
   //                      }
   //                  }
   //                  else
   //                  {
   //                    data[i].push_back
   //                      (std::make_pair(DofObject::invalid_id, Real(0)));
   //                  }
   //              }
   //          };

   //          auto rhss_gather_functor =
   //          [& dof_map,
   //           &_dof_constraints,
   //           &_primal_constraint_values,
   //           max_qoi_num]
   //          (processor_id_type,
   //           const std::vector<dof_id_type> & ids,
   //           std::vector<rhss_datum> & data)
   //          {
   //            // Fill those requests
   //              const std::size_t query_size = ids.size();

   //              data.resize(query_size);
   //              for (std::size_t i=0; i != query_size; ++i)
   //              {
   //                  dof_id_type constrained = ids[i];
   //                  data[i].clear();
   //                  if (_dof_constraints.count(constrained))
   //                  {
   //                      DofConstraintValueMap::const_iterator rhsit =
   //                      _primal_constraint_values.find(constrained);
   //                      data[i].push_back
   //                      ((rhsit == _primal_constraint_values.end()) ?
   //                       0 : rhsit->second);

   //                    // for (unsigned int q = 0; q != max_qoi_num; ++q)
   //                    //   {
   //                        // AdjointDofConstraintValues::const_iterator adjoint_map_it =
   //                        //   _adjoint_constraint_values.find(q);

   //                        // if (adjoint_map_it == _adjoint_constraint_values.end())
   //                        //   {
   //                        //     data[i].push_back(0);
   //                        //     continue;
   //                        //   }

   //                        // const DofConstraintValueMap & constraint_map =
   //                        //   adjoint_map_it->second;

   //                        // DofConstraintValueMap::const_iterator adj_rhsit =
   //                        //   constraint_map.find(constrained);
   //                        // data[i].push_back
   //                        //   ((adj_rhsit == constraint_map.end()) ?
   //                        //    0 : adj_rhsit->second);
   //                      // }
   //                  }
   //              }
   //          };

   //          auto row_action_functor =
   //          [& dof_map,
   //           & _dof_constraints,
   //           & unexpanded_dofs]
   //          (processor_id_type,
   //           const std::vector<dof_id_type> & ids,
   //           const std::vector<row_datum> & data)
   //          {
              
   //              const std::size_t query_size = ids.size();

   //              for (std::size_t i=0; i != query_size; ++i)
   //              {
   //                const dof_id_type constrained = ids[i];

   //                  // An empty row is an constraint with an empty row; for
   //                  // no constraint we use a "no row" placeholder
   //                  if (data[i].empty())
   //                  {
   //                    DofConstraintRow & row = _dof_constraints[constrained];
   //                    row.clear();
   //                  }
   //                  else if (data[i][0].first != DofObject::invalid_id)
   //                  {
   //                      DofConstraintRow & row = _dof_constraints[constrained];
   //                      row.clear();
   //                      for (auto & pair : data[i])
   //                      {
   //                        libmesh_assert_less(pair.first, dof_map.n_dofs());
   //                        row[pair.first] = pair.second;
   //                      }

   //                      // And prepare to check for more recursive constraints
   //                      unexpanded_dofs.insert(constrained);
   //                  }
   //              }
   //          };

   //          auto rhss_action_functor =
   //          [& dof_map,
   //           & _dof_constraints,
   //           & _primal_constraint_values,
   //           max_qoi_num]
   //          (processor_id_type,
   //          const std::vector<dof_id_type> & ids,
   //          const std::vector<rhss_datum> & data)
   //          {
   //              // Add rhs data for any new constraint rows we've found
   //              const std::size_t query_size = ids.size();

   //              for (std::size_t i=0; i != query_size; ++i)
   //              {
   //                if (!data[i].empty())
   //                  {
   //                       dof_id_type constrained = ids[i];
   //                       if (data[i][0] != libMesh::Number(0))
   //                      _primal_constraint_values[constrained] = data[i][0];
   //                       else
   //                      _primal_constraint_values.erase(constrained);

   //                      for (unsigned int q = 0; q != max_qoi_num; ++q)
   //                      {
   //                        // AdjointDofConstraintValues::iterator adjoint_map_it =
   //                        //   _adjoint_constraint_values.find(q);

   //                        // if ((adjoint_map_it == _adjoint_constraint_values.end()) &&
   //                        //     data[i][q+1] == Number(0))
   //                        //   continue;

   //                        // if (adjoint_map_it == _adjoint_constraint_values.end())
   //                        //   adjoint_map_it = _adjoint_constraint_values.insert
   //                        //     (std::make_pair(q,DofConstraintValueMap())).first;

   //                        // DofConstraintValueMap & constraint_map =
   //                        //   adjoint_map_it->second;

   //                        // if (data[i][q+1] != Number(0))
   //                        //   constraint_map[constrained] =
   //                        //     data[i][q+1];
   //                        // else
   //                        //   constraint_map.erase(constrained);
   //                      }
   //                  }
   //              }

   //          };

   //          // Now request constraint rows from other processors
   //          row_datum * row_ex = nullptr;
   //          Parallel::pull_parallel_vector_data
   //          (dof_map.comm(), requested_dof_ids, row_gather_functor,
   //          row_action_functor, row_ex);

   //          //And request constraint right hand sides from other procesors
   //          rhss_datum * rhs_ex = nullptr;
   //          Parallel::pull_parallel_vector_data
   //          (dof_map.comm(), requested_dof_ids, rhss_gather_functor,
   //          rhss_action_functor, rhs_ex);

   //          // We have to keep recursing while the unexpanded set is
   //          // nonempty on *any* processor
   //          unexpanded_set_nonempty = !unexpanded_dofs.empty();
   //          dof_map.comm().max(unexpanded_set_nonempty);

   //      }
   //  }


   //  void Adaptivity::add_constraints_to_send_list(libMesh::DofMap &dof_map, 
   //                                                libMesh::DofConstraints &_dof_constraints,
   //                                                std::vector<libMesh::dof_id_type> &_send_list)
   //  {
  
   //      //std::cout<<"Adaptivity::add_constraints_to_send_list::BEGIN "<<std::endl; 

   //      // Recreate any user or internal constraints
   //      //dof_map.reinit_constraints();


   //      // And clean up the send_list before we first use it
   //      //dof_map.prepare_send_list();


   //      if (dof_map.n_processors() == 1) return;          
        
   //      unsigned int has_constraints = !_dof_constraints.empty();

   //      dof_map.comm().max(has_constraints);

   //      if (!has_constraints) return;

   //      for (const auto & i : _dof_constraints)
   //      {
   //          libMesh::dof_id_type constrained_dof = i.first;

   //          if (!dof_map.local_index(constrained_dof))
   //              continue;

   //          const libMesh::DofConstraintRow & constraint_row = i.second;
              
   //          for (const auto & j : constraint_row)
   //          {
   //              libMesh::dof_id_type constraint_dependency = j.first;

   //              if (dof_map.local_index(constraint_dependency)) continue;

   //              _send_list.push_back(constraint_dependency);
   //          }
   //      }       

   //      //std::cout<<"Adaptivity::add_constraints_to_send_list::END "<<std::endl;  


   //  }


   //  void Adaptivity::scatter_constraints(libMesh::MeshBase & mesh, 
   //                                       libMesh::DofMap &dof_map, 
   //                                       libMesh::DofConstraints &_dof_constraints)
   //  {
           
   //      std::cout<<"Adaptivity::scatter_constraints::Begin "<<std::endl;      
   //      using namespace libMesh;


   //      auto _primal_constraint_values = dof_map.get_primal_constraint_values();

    
   //      // Return immediately if there's nothing to gather
   //      if (dof_map.n_processors() == 1)
   //      return;

   //      // We might get to return immediately if none of the processors
   //      // found any constraints
   //      unsigned int has_constraints = !_dof_constraints.empty()
   //      #ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
   //      || !_node_constraints.empty()
   //      #endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
   //      ;
   //      dof_map.comm().max(has_constraints);
   //      if (!has_constraints)
   //      return;

   //      // We may be receiving packed_range sends out of order with
   //      // parallel_sync tags, so make sure they're received correctly.
   //      Parallel::MessageTag range_tag = dof_map.comm().get_unique_tag();


   //      std::map<processor_id_type, std::set<dof_id_type>> pushed_ids;

   //      // Collect the dof constraints I need to push to each processor
   //      dof_id_type constrained_proc_id = 0;

   //      for (auto & i : _dof_constraints)
   //      {
   //         const dof_id_type constrained = i.first;
           
   //         while (constrained >= dof_map.end_dof(constrained_proc_id))
   //          constrained_proc_id++;

   //          if (constrained_proc_id != dof_map.processor_id())
   //              continue;

   //          DofConstraintRow & row = i.second;
   //          for (auto & j : row)
   //          {
   //            const dof_id_type constraining = j.first;

   //            processor_id_type constraining_proc_id = 0;
   //            while (constraining >= dof_map.end_dof(constraining_proc_id))
   //              constraining_proc_id++;

   //            if (constraining_proc_id != dof_map.processor_id() &&
   //                constraining_proc_id != constrained_proc_id)
   //              pushed_ids[constraining_proc_id].insert(constrained);
   //          }
   //      }

   //      // Pack the dof constraint rows and rhs's to push

   //      std::map<processor_id_type,
   //            std::vector<std::vector<std::pair<dof_id_type, Real>>>>
   //      pushed_keys_vals, pushed_keys_vals_to_me;

   //      std::map<processor_id_type, std::vector<std::pair<dof_id_type, libMesh::Number>>>
   //      pushed_ids_rhss, pushed_ids_rhss_to_me;

   //      auto gather_ids =
   //      [& dof_map,
   //       & pushed_ids,
   //       & _primal_constraint_values,
   //       & _dof_constraints,
   //       & pushed_keys_vals,
   //       & pushed_ids_rhss]
   //      ()
   //      {
   //        for (auto & pid_id_pair : pushed_ids)
   //          {
   //            const processor_id_type pid = pid_id_pair.first;
   //            const std::set<dof_id_type>
   //              & pid_ids = pid_id_pair.second;

   //            const std::size_t ids_size = pid_ids.size();
   //            std::vector<std::vector<std::pair<dof_id_type, Real>>> &
   //              keys_vals = pushed_keys_vals[pid];
   //            std::vector<std::pair<dof_id_type,libMesh::Number>> &
   //              ids_rhss = pushed_ids_rhss[pid];
   //            keys_vals.resize(ids_size);
   //            ids_rhss.resize(ids_size);

   //            std::size_t push_i;
   //            std::set<dof_id_type>::const_iterator it;
   //            for (push_i = 0, it = pid_ids.begin();
   //                 it != pid_ids.end(); ++push_i, ++it)
   //              {
   //                const dof_id_type constrained = *it;
   //                DofConstraintRow & row = _dof_constraints[constrained];
   //                keys_vals[push_i].assign(row.begin(), row.end());

   //                DofConstraintValueMap::const_iterator rhsit =
   //                  _primal_constraint_values.find(constrained);
   //                ids_rhss[push_i].first = constrained;
   //                ids_rhss[push_i].second =
   //                  (rhsit == _primal_constraint_values.end()) ?
   //                  0 : rhsit->second;
   //              }
   //          }
   //      };

   //      gather_ids();

   //      auto ids_rhss_action_functor =
   //      [& pushed_ids_rhss_to_me]
   //      (processor_id_type pid,
   //       const std::vector<std::pair<dof_id_type, libMesh::Number>> & data)
   //      {
   //        pushed_ids_rhss_to_me[pid] = data;
   //      };

   //      auto keys_vals_action_functor =
   //      [& pushed_keys_vals_to_me]
   //      (processor_id_type pid,
   //       const std::vector<std::vector<std::pair<dof_id_type, Real>>> & data)
   //      {
   //        pushed_keys_vals_to_me[pid] = data;
   //      };

   //      Parallel::push_parallel_vector_data
   //      (dof_map.comm(), pushed_ids_rhss, ids_rhss_action_functor);
   //      Parallel::push_parallel_vector_data
   //      (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

   //      // Now work on traded dof constraint rows
   //      auto receive_dof_constraints =
   //      [& dof_map,
   //       & _primal_constraint_values,
   //       & _dof_constraints,
   //       & pushed_ids_rhss_to_me,
   //       & pushed_keys_vals_to_me]
   //      ()
   //      {
   //        for (auto & pid_id_pair : pushed_ids_rhss_to_me)
   //          {
   //            const processor_id_type pid = pid_id_pair.first;
   //            const auto & ids_rhss = pid_id_pair.second;
   //            const auto & keys_vals = pushed_keys_vals_to_me[pid];

   //            libmesh_assert_equal_to
   //              (ids_rhss.size(), keys_vals.size());

   //            // Add the dof constraints that I've been sent
   //            for (auto i : index_range(ids_rhss))
   //              {
   //                dof_id_type constrained = ids_rhss[i].first;

   //                // If we don't already have a constraint for & dof_map dof,
   //                // add the one we were sent
   //                if (!dof_map.is_constrained_dof(constrained))
   //                  {
   //                    DofConstraintRow & row = _dof_constraints[constrained];
   //                    for (auto & key_val : keys_vals[i])
   //                      {
   //                        libmesh_assert_less(key_val.first, dof_map.n_dofs());
   //                        row[key_val.first] = key_val.second;
   //                      }
   //                    if (ids_rhss[i].second != libMesh::Number(0))
   //                      _primal_constraint_values[constrained] =
   //                        ids_rhss[i].second;
   //                    else
   //                      _primal_constraint_values.erase(constrained);
   //                  }
   //              }
   //          }
   //      };

   //      receive_dof_constraints();


   //      typedef std::map<dof_id_type, std::set<dof_id_type>> DofConstrainsMap;
   //      DofConstrainsMap dof_id_constrains;

   //      for (auto & i : _dof_constraints)
   //      {
   //        const dof_id_type constrained = i.first;
   //        DofConstraintRow & row = i.second;
   //        for (const auto & j : row)
   //          {
   //            const dof_id_type constraining = j.first;

   //            dof_id_type constraining_proc_id = 0;
   //            while (constraining >= dof_map.end_dof(constraining_proc_id))
   //              constraining_proc_id++;

   //            if (constraining_proc_id == dof_map.processor_id())
   //              dof_id_constrains[constraining].insert(constrained);
   //          }
   //      }

   //      // Loop over all foreign elements, find any supporting our
   //      // constrained dof indices.
   //      pushed_ids.clear();

   //      for (const auto & elem : as_range(mesh.active_not_local_elements_begin(),
   //                                      mesh.active_not_local_elements_end()))
   //      {
   //        std::vector<dof_id_type> my_dof_indices;
   //        dof_map.dof_indices (elem, my_dof_indices);

   //        for (const auto & dof : my_dof_indices)
   //          {
   //            DofConstrainsMap::const_iterator dcmi = dof_id_constrains.find(dof);
   //            if (dcmi != dof_id_constrains.end())
   //              {
   //                for (const auto & constrained : dcmi->second)
   //                  {
   //                    dof_id_type the_constrained_proc_id = 0;
   //                    while (constrained >= dof_map.end_dof(the_constrained_proc_id))
   //                      the_constrained_proc_id++;

   //                    const processor_id_type elemproc = elem->processor_id();
   //                    if (elemproc != the_constrained_proc_id)
   //                      pushed_ids[elemproc].insert(constrained);
   //                  }
   //              }
   //          }
   //      }

   //      pushed_ids_rhss.clear();
   //      pushed_ids_rhss_to_me.clear();
   //      pushed_keys_vals.clear();
   //      pushed_keys_vals_to_me.clear();

   //      gather_ids();

   //      // Trade pushed dof constraint rows
   //      Parallel::push_parallel_vector_data
   //      (dof_map.comm(), pushed_ids_rhss, ids_rhss_action_functor);
   //      Parallel::push_parallel_vector_data
   //      (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

   //      receive_dof_constraints();

   //      // Finally, we need to handle the case of remote dof coupling.  If a
   //      // processor's element is coupled to a ghost element, then the
   //      // processor needs to know about all constraints which affect the
   //      // dofs on that ghost element, so we'll have to query the ghost
   //      // element's owner.

   //      GhostingFunctor::map_type elements_to_couple;

   //      // Man, I wish we had guaranteed unique_ptr availability...
   //      std::set<CouplingMatrix*> temporary_coupling_matrices;

   //      Adaptivity::merge_ghost_functor_outputs(elements_to_couple,temporary_coupling_matrices,
   //                                              dof_map.coupling_functors_begin(),
   //                                              dof_map.coupling_functors_end(),
   //                                              mesh.active_local_elements_begin(),
   //                                              mesh.active_local_elements_end(),
   //                                              dof_map.processor_id());

   //      std::set<dof_id_type> requested_dofs;

   //      for (const auto & pr : elements_to_couple)
   //      {
   //        const Elem * elem = pr.first;

   //        std::vector<dof_id_type> element_dofs;

   //        dof_map.dof_indices(elem, element_dofs);

   //        for (auto dof : element_dofs)requested_dofs.insert(dof);
   //      } 

   //      gather_constraints(mesh, requested_dofs, _dof_constraints, dof_map, false);
   // }  


   //  void Adaptivity::merge_ghost_functor_outputs(libMesh::GhostingFunctor::map_type & elements_to_ghost,
   //                          std::set<libMesh::CouplingMatrix *> & temporary_coupling_matrices,
   //                          const std::set<libMesh::GhostingFunctor *>::iterator & gf_begin,
   //                          const std::set<libMesh::GhostingFunctor *>::iterator & gf_end,
   //                          const libMesh::MeshBase::const_element_iterator & elems_begin,
   //                          const libMesh::MeshBase::const_element_iterator & elems_end,
   //                          libMesh::processor_id_type p)
   // {
   //      using namespace libMesh;

   //      for (const auto & gf : as_range(gf_begin, gf_end))
   //      {
   //        GhostingFunctor::map_type more_elements_to_ghost;

   //        libmesh_assert(gf);
   //        (*gf)(elems_begin, elems_end, p, more_elements_to_ghost);

   //          for (const auto & pr : more_elements_to_ghost)
   //          {
   //            GhostingFunctor::map_type::iterator existing_it =
   //              elements_to_ghost.find (pr.first);
   //              if (existing_it == elements_to_ghost.end())
   //              elements_to_ghost.insert(pr);
   //              else
   //              {
   //                  if (existing_it->second)
   //                  {
   //                      if (pr.second)
   //                      {
                         
   //                          if (temporary_coupling_matrices.empty() ||
   //                             temporary_coupling_matrices.find(const_cast<CouplingMatrix *>(existing_it->second)) == temporary_coupling_matrices.end())
   //                          {
   //                            CouplingMatrix * cm = new CouplingMatrix(*existing_it->second);
   //                            temporary_coupling_matrices.insert(cm);
   //                            existing_it->second = cm;
   //                          }
   //                          const_cast<CouplingMatrix &>(*existing_it->second) &= *pr.second;
   //                      }
   //                      else
   //                      {
   //                          std::set<CouplingMatrix *>::iterator temp_it =
   //                          temporary_coupling_matrices.find(const_cast<CouplingMatrix *>(existing_it->second));
   //                          if (temp_it != temporary_coupling_matrices.end())
   //                             temporary_coupling_matrices.erase(temp_it);

   //                          existing_it->second = nullptr;
   //                      }
   //                  }               
   //              }
   //          }
   //      }
   //  }


   //  void Adaptivity::check_for_constraint_loops(libMesh::DofMap &dof_map, libMesh::DofConstraints &_dof_constraints)
   //  {
        

   //      using namespace libMesh;

   //      typedef std::set<dof_id_type> RCSet;
   //      RCSet unexpanded_set;

   //      DofConstraints dof_constraints_copy = _dof_constraints;

   //      for (const auto & i : dof_constraints_copy)
   //      unexpanded_set.insert(i.first);

   //      while (!unexpanded_set.empty())
   //      for (RCSet::iterator i = unexpanded_set.begin();
   //           i != unexpanded_set.end(); /* nothing */)
   //      {
   //          // If the DOF is constrained
   //          DofConstraints::iterator
   //          pos = dof_constraints_copy.find(*i);

   //          libmesh_assert (pos != dof_constraints_copy.end());

   //          DofConstraintRow & constraint_row = pos->second;

   //          // Comment out "rhs" parts of & dof_map method copied from process_constraints
   //          // DofConstraintValueMap::iterator rhsit =
   //          //   _primal_constraint_values.find(*i);
   //          // Number constraint_rhs = (rhsit == _primal_constraint_values.end()) ?
   //          //   0 : rhsit->second;

   //          std::vector<dof_id_type> constraints_to_expand;

   //          for (const auto & item : constraint_row)
   //          if (item.first != *i && dof_map.is_constrained_dof(item.first))
   //          {
   //            unexpanded_set.insert(item.first);
   //            constraints_to_expand.push_back(item.first);
   //          }

   //          for (const auto & expandable : constraints_to_expand)
   //          {
   //          const Real this_coef = constraint_row[expandable];

   //          DofConstraints::const_iterator
   //            subpos = dof_constraints_copy.find(expandable);

   //          libmesh_assert (subpos != dof_constraints_copy.end());

   //          const DofConstraintRow & subconstraint_row = subpos->second;

   //          for (const auto & item : subconstraint_row)
   //            {
   //              if (item.first == expandable)
   //                libmesh_error_msg("Constraint loop detected");

   //              constraint_row[item.first] += item.second * this_coef;
   //            }


   //          constraint_row.erase(expandable);
   //      }

   //          // Comment out "rhs" parts of & dof_map method copied from process_constraints
   //          // if (rhsit == _primal_constraint_values.end())
   //          //   {
   //          //     if (constraint_rhs != Number(0))
   //          //       _primal_constraint_values[*i] = constraint_rhs;
   //          //     else
   //          //       _primal_constraint_values.erase(*i);
   //          //   }
   //          // else
   //          //   {
   //          //     if (constraint_rhs != Number(0))
   //          //       rhsit->second = constraint_rhs;
   //          //     else
   //          //       _primal_constraint_values.erase(rhsit);
   //          //   }

   //          if (constraints_to_expand.empty())
   //            i = unexpanded_set.erase(i);
   //          else
   //            ++i;
   //      }
   //  }


          
}