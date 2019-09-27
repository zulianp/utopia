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



namespace utopia {

    void Adaptivity::compute_all_constraints(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        libMesh::DofConstraints &constraints)
    {

      
        using uint = unsigned int;
        
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

                    compute_constraints(constraints, dof_map, var_num, elem, mesh.mesh_dimension());
                }
            }
        }

        libMesh::MeshBase &mesh_copy=const_cast<libMesh::MeshBase&>(mesh);

        libMesh::DofMap &dof_copy=const_cast<libMesh::DofMap&>(dof_map);

        uint n_vars = dof_map.n_variables();

        for(uint var_num = 0; var_num < n_vars; ++var_num) 
        {

            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            fe_type.order = static_cast<libMesh::Order>(fe_type.order);

            if (fe_type.order>0)
            {
                process_constraints(mesh_copy, dof_copy, constraints);
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

        dof_constraints_.clear();
        

        libMesh::FEType fe_type = dof_map.variable_type(var_num);

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
                    compute_constraints(dof_constraints_, dof_map,  var_num, elem, mesh.mesh_dimension());
                }
                
            }

            libMesh::MeshBase &mesh_copy=const_cast<libMesh::MeshBase&>(mesh);

            libMesh::DofMap &dof_copy=const_cast<libMesh::DofMap&>(dof_map);

            Adaptivity::process_constraints(mesh_copy, dof_copy, dof_constraints_);
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
        constraint_matrix(mesh, dof_map, var_num, M, S);

    }

    void Adaptivity::constraint_matrix(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, int var_num, USparseMatrix &M, USparseMatrix &S)
    {
        using IndexArray = Traits<USparseMatrix>::IndexArray;
        
        std::vector<SizeType> index;

        assemble_constraint(mesh, dof_map, var_num);
        
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
                    
                    libmesh_assert (pos != dof_constraints_.end());
                    
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

    //   std::cout<<"Adaptivity::process_constraints::BEGIN "<<std::endl;

       
       std::vector<int> index; 

       compute_boundary_nodes(mesh, dof_map, 0,0, index);

       auto on_boundary = libMesh::MeshTools::find_boundary_nodes(mesh);

       typedef std::set<libMesh::dof_id_type> RCSet;

       RCSet unexpanded_set;

       for (const auto & i : _dof_constraints)
       { 
           unexpanded_set.insert(i.first);
       }

       while (!unexpanded_set.empty())
        for (RCSet::iterator i = unexpanded_set.begin(); i != unexpanded_set.end(); /* nothing */)
        {
                              
            libMesh::DofConstraints::iterator pos = _dof_constraints.find(*i);

            libmesh_assert( pos != _dof_constraints.end());

            libMesh::DofConstraintRow & constraint_row = pos->second;

            std::vector<libMesh::dof_id_type> constraints_to_expand;

            for (const auto & item : constraint_row)
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
         
                    if (check ==true) {
                          unexpanded_set.insert(item.first);     

                          constraints_to_expand.push_back(item.first);
                    }
                }

            for (const auto & expandable : constraints_to_expand)
            {
                const double this_coef = constraint_row[expandable];

                libMesh::DofConstraints::const_iterator subpos = _dof_constraints.find(expandable);

                libmesh_assert(subpos != _dof_constraints.end());

                const libMesh::DofConstraintRow & subconstraint_row = subpos->second;

                for (const auto & item : subconstraint_row)
                  {
                    // Assert that the constraint does not form a cycle.
                    libmesh_assert(item.first != expandable);
                    
                    constraint_row[item.first] += item.second * this_coef;
                    
                    //std::cout<<"item.second * this_coef=>"<<item.second * this_coef<<std::endl;
                  }

                constraint_row.erase(expandable);
            }

            if (constraints_to_expand.empty()) 
            {
                i = unexpanded_set.erase(i);
            }
            else
            {
              ++i;
            }
        }
        
    //  std::cout<<"Adaptivity::process_constraints::END "<<std::endl;    

    }

    void Adaptivity::allgather_recursive_constraints(libMesh::MeshBase & mesh, 
                                                    libMesh::DofConstraints &_dof_constraints, 
                                                    libMesh::DofMap &dof_map)
    {
       // std::cout<<"Adaptivity::allgather_recursive_constraints::BEGIN "<<std::endl;


      if (dof_map.n_processors() == 1)
        return;


      unsigned int has_constraints = !_dof_constraints.empty();

      dof_map.comm().max(has_constraints);

      if (!has_constraints)
        return;
      {
        std::map<libMesh::processor_id_type, std::set<libMesh::dof_id_type>> pushed_ids;



        const unsigned int sys_num = dof_map.sys_number();

        // Collect the constraints to push to each processor
        for (auto & elem : as_range(mesh.active_not_local_elements_begin(),
                                    mesh.active_not_local_elements_end()))
          {
            const unsigned short n_nodes = elem->n_nodes();

            {
              const unsigned int n_vars = elem->n_vars(sys_num);
              for (unsigned int v=0; v != n_vars; ++v)
                {
                  const unsigned int n_comp = elem->n_comp(sys_num,v);
                  for (unsigned int c=0; c != n_comp; ++c)
                    {
                      const unsigned int id =
                        elem->dof_number(sys_num,v,c);
                      if (dof_map.is_constrained_dof(id))
                        pushed_ids[elem->processor_id()].insert(id);
                    }
                }
            }

            for (unsigned short n = 0; n != n_nodes; ++n)
              {
                const libMesh::Node & node = elem->node_ref(n);
                const unsigned int n_vars = node.n_vars(sys_num);
                for (unsigned int v=0; v != n_vars; ++v)
                  {
                    const unsigned int n_comp = node.n_comp(sys_num,v);
                    for (unsigned int c=0; c != n_comp; ++c)
                      {
                        const unsigned int id =
                          node.dof_number(sys_num,v,c);
                        if (dof_map.is_constrained_dof(id))
                          pushed_ids[elem->processor_id()].insert(id);
                      }
                  }
              }

          }

        // Rewrite those id sets as vectors for sending and receiving,
        // then find the corresponding data for each id, then push it all.
        std::map<libMesh::processor_id_type, std::vector<libMesh::dof_id_type>>
          pushed_id_vecs, received_id_vecs;
        for (auto & p : pushed_ids)
          pushed_id_vecs[p.first].assign(p.second.begin(), p.second.end());

        std::map<libMesh::processor_id_type, std::vector<std::vector<std::pair<libMesh::dof_id_type,double>>>>
          pushed_keys_vals, received_keys_vals;

        for (auto & p : pushed_id_vecs)
        {
            auto & keys_vals = pushed_keys_vals[p.first];
            keys_vals.reserve(p.second.size());

            for (auto & pushed_id : p.second)
              {
                const libMesh::DofConstraintRow & row = _dof_constraints[pushed_id];
                keys_vals.emplace_back(row.begin(), row.end());

              }
        }

        auto ids_action_functor =
          [& received_id_vecs]
          (libMesh::processor_id_type pid,
           const std::vector<libMesh::dof_id_type> & data)
          {
            received_id_vecs[pid] = data;
          };

        libMesh::Parallel::push_parallel_vector_data
          (dof_map.comm(), pushed_id_vecs, ids_action_functor);

        auto keys_vals_action_functor =
          [& received_keys_vals]
          (libMesh::processor_id_type pid,
           const std::vector<std::vector<std::pair<libMesh::dof_id_type,double>>> & data)
          {
            received_keys_vals[pid] = data;
          };

        libMesh::Parallel::push_parallel_vector_data
          (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);



        for (auto & p : received_id_vecs)
          {
            const libMesh::processor_id_type pid = p.first;
            const auto & pushed_ids_to_me = p.second;
            libmesh_assert(received_keys_vals.count(pid));

            const auto & pushed_keys_vals_to_me = received_keys_vals.at(pid);


            libmesh_assert_equal_to (pushed_ids_to_me.size(),
                                     pushed_keys_vals_to_me.size());


            for (auto i : libMesh::index_range(pushed_ids_to_me))
              {
                libMesh::dof_id_type constrained = pushed_ids_to_me[i];

        
                if (!dof_map.is_constrained_dof(constrained))
                  {
                    libMesh::DofConstraintRow & row = _dof_constraints[constrained];
                    for (auto & kv : pushed_keys_vals_to_me[i])
                      {
                        libmesh_assert_less(kv.first, dof_map.n_dofs());
                        row[kv.first] = kv.second;
                      }
                  }
              }
          }
      }


      typedef std::set<libMesh::dof_id_type> DoF_RCSet;

      DoF_RCSet unexpanded_dofs;

      for (const auto & i : _dof_constraints) unexpanded_dofs.insert(i.first);


      // std::cout<<"Adaptivity::allgather_recursive_constraints::END"<<std::endl;


      gather_constraints(mesh, unexpanded_dofs, _dof_constraints, dof_map, false);

      
    }



    void Adaptivity::gather_constraints (libMesh::MeshBase & mesh,
                                         std::set<libMesh::dof_id_type> & unexpanded_dofs, 
                                         libMesh::DofConstraints & _dof_constraints,
                                         libMesh::DofMap &dof_map,
                                         bool look_for_constrainees)
    {
      
      //  std::cout<<"Adaptivity::gather_constraints::BEGIN "<<std::endl;  


        typedef std::set<libMesh::dof_id_type> DoF_RCSet;

        bool unexpanded_set_nonempty = !unexpanded_dofs.empty();

        dof_map.comm().max(unexpanded_set_nonempty);

        while (unexpanded_set_nonempty)
        {
              
              DoF_RCSet   dof_request_set;


              std::map<libMesh::processor_id_type, std::vector<libMesh::dof_id_type>> requested_dof_ids;


              std::map<libMesh::processor_id_type, libMesh::dof_id_type> dof_ids_on_proc;


            for (const auto & unexpanded_dof : unexpanded_dofs)
            {
                  libMesh::DofConstraints::const_iterator
                    pos = _dof_constraints.find(unexpanded_dof);

       
                if (pos == _dof_constraints.end())
                {
                      if (!dof_map.local_index(unexpanded_dof) &&
                          !_dof_constraints.count(unexpanded_dof) )
                        dof_request_set.insert(unexpanded_dof);
                }

                else
                {
                      const libMesh::DofConstraintRow & row = pos->second;
                     
                    for (const auto & j : row)
                    {
                          const libMesh::dof_id_type constraining_dof = j.first;


                          if (!dof_map.local_index(constraining_dof) &&
                              !_dof_constraints.count(constraining_dof))
                            dof_request_set.insert(constraining_dof);
                    }
                }
            }

              unexpanded_dofs.clear();


              libMesh::processor_id_type proc_id = 0;
              for (const auto & i : dof_request_set)
              {
                  while (i >= dof_map.end_dof(proc_id))
                  {
                    proc_id++;
                    dof_ids_on_proc[proc_id]++;
                  }
              }

              for (auto & pair : dof_ids_on_proc)
              {
                  requested_dof_ids[pair.first].reserve(pair.second);
              }

           
              proc_id = 0;

              for (const auto & i : dof_request_set)
              {
                  while (i >= dof_map.end_dof(proc_id))
                  {
                    proc_id++;
                    requested_dof_ids[proc_id].push_back(i);
                }
              }

              unexpanded_set_nonempty = !unexpanded_dofs.empty();
              dof_map.comm().max(unexpanded_set_nonempty);
            }



     //   std::cout<<"Adaptivity::gather_constraints::END "<<std::endl; 


    }

    void Adaptivity::add_constraints_to_send_list(libMesh::DofMap &dof_map, 
                                                  libMesh::DofConstraints &_dof_constraints)
    {
      
        //std::cout<<"Adaptivity::add_constraints_to_send_list::BEGIN "<<std::endl; 


        if (dof_map.n_processors() == 1) return;

          
          unsigned int has_constraints = !_dof_constraints.empty();

          dof_map.comm().max(has_constraints);

        if (!has_constraints) return;

        for (const auto & i : _dof_constraints)
        {
              libMesh::dof_id_type constrained_dof = i.first;

              if (!dof_map.local_index(constrained_dof))
                continue;

              const libMesh::DofConstraintRow & constraint_row = i.second;
              
            for (const auto & j : constraint_row)
            {
                  libMesh::dof_id_type constraint_dependency = j.first;

                  if (dof_map.local_index(constraint_dependency)) continue;

                   std::vector<libMesh::dof_id_type> _send_list = dof_map.get_send_list();

                  _send_list.push_back(constraint_dependency);
            }
        }       

        //std::cout<<"Adaptivity::add_constraints_to_send_list::END "<<std::endl;  


    }


    void Adaptivity::scatter_constraints(libMesh::MeshBase & mesh, 
                                         libMesh::DofMap &dof_map, 
                                         libMesh::DofConstraints &_dof_constraints)
    {




      // std::cout<<"Adaptivity::scatter_constraints::BEGIN"<<std::endl;

      // This function must be run on all processors at once
      //parallel_object_only();

      // Return immediately if there's nothing to gather
      if (dof_map.n_processors() == 1)
        return;

      // We might get to return immediately if none of the processors
      // found any constraints
      unsigned int has_constraints = !_dof_constraints.empty();

        
      dof_map.comm().max(has_constraints);

      if (!has_constraints)
        return;

      // libMesh::Parallel::MessageTag range_tag = dof_map.comm().get_unique_tag();


      std::map<libMesh::processor_id_type, std::set<libMesh::dof_id_type>> pushed_ids;

      // Collect the dof constraints I need to push to each processor
      libMesh::dof_id_type constrained_proc_id = 0;

    for (auto & i : _dof_constraints)
    {
          const libMesh::dof_id_type constrained = i.first;
          while (constrained >= dof_map.end_dof(constrained_proc_id))
            constrained_proc_id++;

          if (constrained_proc_id != dof_map.processor_id())
            continue;

          libMesh::DofConstraintRow & row = i.second;
          for (auto & j : row)
            {
              const libMesh::dof_id_type constraining = j.first;

              libMesh::processor_id_type constraining_proc_id = 0;
              while (constraining >= dof_map.end_dof(constraining_proc_id))
                constraining_proc_id++;

              if (constraining_proc_id != dof_map.processor_id() &&
                  constraining_proc_id != constrained_proc_id)
                pushed_ids[constraining_proc_id].insert(constrained);
            }
        }

      // Pack the dof constraint rows and rhs's to push

      std::map<libMesh::processor_id_type,
              std::vector<std::vector<std::pair<libMesh::dof_id_type, double>>>>
        pushed_keys_vals, pushed_keys_vals_to_me;


      auto keys_vals_action_functor =
        [& pushed_keys_vals_to_me]
        (libMesh::processor_id_type pid,
         const std::vector<std::vector<std::pair<libMesh::dof_id_type, double>>> & data)
        {
          pushed_keys_vals_to_me[pid] = data;
        };


      libMesh::Parallel::push_parallel_vector_data
        (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

      

      typedef std::map<libMesh::dof_id_type, std::set<libMesh::dof_id_type>> DofConstrainsMap;
      DofConstrainsMap dof_id_constrains;

      for (auto & i : _dof_constraints)
        {
          const libMesh::dof_id_type constrained = i.first;
          libMesh::DofConstraintRow & row = i.second;
          for (const auto & j : row)
            {
              const libMesh::dof_id_type constraining = j.first;

              libMesh::dof_id_type constraining_proc_id = 0;
              while (constraining >= dof_map.end_dof(constraining_proc_id))
                constraining_proc_id++;

              if (constraining_proc_id == dof_map.processor_id())
                dof_id_constrains[constraining].insert(constrained);
            }
        }

      // Loop over all foreign elements, find any supporting our
      // constrained dof indices.
      pushed_ids.clear();

      for (const auto & elem : as_range(mesh.active_not_local_elements_begin(),
                                        mesh.active_not_local_elements_end()))
        {
          std::vector<libMesh::dof_id_type> my_dof_indices;
          dof_map.dof_indices (elem, my_dof_indices);

          for (const auto & dof : my_dof_indices)
            {
              DofConstrainsMap::const_iterator dcmi = dof_id_constrains.find(dof);
              if (dcmi != dof_id_constrains.end())
                {
                  for (const auto & constrained : dcmi->second)
                    {
                      libMesh::dof_id_type the_constrained_proc_id = 0;
                      while (constrained >= dof_map.end_dof(the_constrained_proc_id))
                        the_constrained_proc_id++;

                      const libMesh::processor_id_type elemproc = elem->processor_id();
                      if (elemproc != the_constrained_proc_id)
                        pushed_ids[elemproc].insert(constrained);
                    }
                }
            }
        }


      pushed_keys_vals.clear();
      pushed_keys_vals_to_me.clear();

      libMesh::Parallel::push_parallel_vector_data
        (dof_map.comm(), pushed_keys_vals, keys_vals_action_functor);

    }


    

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

                            index.push_back(node_dof);
                        }
                    }
                }
            }
        }


    // std::cout<<"Adaptivity::compute_boundary_nodes::END "<<std::endl; 
    }
        
}