#include "utopia_LibMeshBackend.hpp"
#include "libmesh/petsc_vector.h"
#include "utopia_Adaptivity.hpp"
#include "libmesh/remote_elem.h"
#include "libmesh/fe_interface.h"
namespace utopia {

    void apply_boundary_conditions(LibMeshFunctionSpace &V,
                                  USparseMatrix &mat, UVector &vec)
    {
        using SizeType = Traits<UVector>::SizeType;

        if(utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions Adaptivity begin: "  << std::endl;
        }

        Chrono c;
        c.start();

        assert(!empty(mat));
        assert(!empty(vec));

       std::vector<SizeType> index, dirichel_id;
       std::vector<SizeType> index_local;

       auto on_boundary = libMesh::MeshTools::find_boundary_nodes(V.mesh());      

       auto & mesh = V.mesh();

       auto & dof_map = V.dof_map();

       libMesh::DofConstraintValueMap &rhs_values = V.dof_map().get_primal_constraint_values();  


      {
          libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();
          
          const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();

          std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;
         
          std::unique_ptr<const libMesh::Elem> my_side, parent_side_0;

          libMesh::FEType fe_type = dof_map.variable_type(0);
      

          
          for ( ; it != end_it; ++it)
          {
            const libMesh::Elem * ele = *it; 

            //fe_type.order = static_cast<libMesh::Order>(fe_type.order + ele->p_level());

            const auto *ele_parent_0 = ele->top_parent();

            for(int jj=0; jj<ele_parent_0->n_sides(); jj++) 
            {

             
              
              libmesh_assert(ele_parent_0);

              //auto my_side = ele->build_side_ptr(jj);

              //std::cout<<"my_side"<<*my_side<<'\n'; 



              auto parent_side_0 = ele_parent_0->build_side_ptr(jj);


             

              // my_dof_indices.reserve (my_side->n_nodes()); 
            
              // parent_dof_indices.reserve (parent_side_0->n_nodes());

              // my_dof_indices.clear();

              // parent_dof_indices.clear();

              // V.dof_map().dof_indices (my_side.get(), my_dof_indices,  0);

              // V.dof_map().dof_indices (parent_side_0.get(), parent_dof_indices, 0);

              // const unsigned int n_side_dofs = libMesh::FEInterface::n_dofs(V.mesh().mesh_dimension()-1, fe_type, my_side->type());
                  
              // const unsigned int n_parent_side_dofs = libMesh::FEInterface::n_dofs(V.mesh().mesh_dimension()-1, fe_type, parent_side_0->type());


              index_local.clear();

              for (int ll=0; ll<parent_side_0->n_nodes(); ll++)
              {
            
                const libMesh::Node * node_0 = parent_side_0->node_ptr(ll);

                //auto node_dof_0 = parent_dof_indices[ll];


                const libMesh::dof_id_type node_dof_0 = node_0->dof_number(V.equation_system().number(), 0, 0); 
               
                if(dof_map.is_constrained_dof(node_dof_0)) {
                    
                    index_local.push_back(node_dof_0);

                    auto valpos = rhs_values.find(node_dof_0);

                    index.push_back(node_dof_0);
                    
                                            

                }
              }

               if(index_local.size()==parent_side_0->n_nodes()){

                auto bc_id = mesh.get_boundary_info().boundary_id(ele_parent_0,jj);

                std::cout<<"bc_id"<<bc_id<<std::endl;

                std::cout<<"my_side_p"<<*parent_side_0<<'\n';  

                auto check = (std::find(dirichel_id.begin(), dirichel_id.end(), bc_id) != dirichel_id.end());


                if (!check) dirichel_id.push_back(bc_id);

                //auto side = ele->build_side_ptr(jj);
             }
           }
         }
       }
    
      //  std::cout<<"dirichel_id"<<dirichel_id.size()<<std::endl;

       {
            libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();
            
            const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();
            
            for ( ; it != end_it; ++it)
            {
                const libMesh::Elem * ele = *it;

                for(int kk=0; kk<ele->n_sides(); kk++) {     

                    auto neigh = ele->neighbor_ptr(kk); 

                    auto bc_id = mesh.get_boundary_info().boundary_id(ele,kk);

                    auto check = (std::find(dirichel_id.begin(), dirichel_id.end(), bc_id) != dirichel_id.end());

                    if(check)
                    {
                         index_local.clear();

                         auto side = ele->build_side_ptr(kk);

                        for (int ll=0; ll<side->n_nodes(); ll++)
                        {
                          
                            const libMesh::Node * node = side->node_ptr(ll);

                            const libMesh::dof_id_type node_dof = node->dof_number(V.equation_system().number(), 0, 0); 

                            index.push_back(node_dof);


                        //     auto valpos = rhs_values.find(node_dof);

                        //     double value = valpos->second;

                        //     std::cout<<"value "<<value<<std::endl;
                   

                        //     if(dof_map.is_constrained_dof(node_dof)) 
                        //     {
                                   
                        //         index_local.push_back(node_dof);
           
                        //     }

                        }

                        // if(index_local.size()==side->n_nodes())
                        // {

                        //    index.insert(index.end(), index_local.begin(), index_local.end());

                        // }
                    }
                }
            }
        }
        

        const bool has_constaints = V.dof_map().constraint_rows_begin() != V.dof_map().constraint_rows_end();


        Size ls = local_size(mat);

        Size s = size(mat);

        set_zero_rows(mat, index, 1.);

        Write<UVector> w_v(vec, utopia::GLOBAL_INSERT);

        std::vector<SizeType> I(1,0);
        std::vector<double> value(1, 0);

        //utopia::disp(vec);

        if(has_constaints) 
        {
            libMesh::DofConstraintValueMap &rhs_values = V.dof_map().get_primal_constraint_values();

            //rhs_values.

            Range r = range(vec);

            for(auto it=index.begin(); it < index.end(); ++it)
            {
              int i = *it;
              auto valpos = rhs_values.find(i);
              I[0] = i;
              value[0]=valpos->second;
              vec.set(I, value);

          }
        }
        
        //utopia::disp(vec);

        c.stop();

        if(utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions end: " << c << std::endl;
        }

        // std::cout << "apply_boundary_conditions end: " << c << std::endl;
    }


     void apply_boundary_conditions(libMesh::DofMap &dof_map, USparseMatrix &mat, UVector &vec)
    {
        if(utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions begin: "  << std::endl;
        }

        Chrono c;
        c.start();

        assert(!empty(mat));
        assert(!empty(vec));

        using SizeType = Traits<UVector>::SizeType;

        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();


        Size ls = local_size(mat);
        Size s = size(mat);

        std::vector<SizeType> index;

        Range rr = range(vec);

        if(has_constaints) {
            for(SizeType i = rr.begin(); i < rr.end(); ++i) {
                if( dof_map.is_constrained_dof(i)) {
                    index.push_back(i);
                }
            }
        }

        set_zero_rows(mat, index, 1.);

        Write<UVector> w_v(vec);

        if(has_constaints) {
            libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

            Range r = range(vec);
            for(SizeType i = r.begin(); i < r.end(); ++i) {
                if(dof_map.is_constrained_dof(i)) {
                    auto valpos = rhs_values.find(i);
                    vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                }
            }
        }

        c.stop();

        if(utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions end: " << c << std::endl;
        }
    }


    void apply_boundary_conditions(libMesh::DofMap &dof_map, UVector &vec)
    {
        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Write<UVector> w_v(vec);

        if(has_constaints) {
            libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

            Range r = range(vec);
            for(SizeType i = r.begin(); i < r.end(); ++i) {
                if(dof_map.is_constrained_dof(i)) {
                    auto valpos = rhs_values.find(i);
                    vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                }
            }
        }

    }

}



