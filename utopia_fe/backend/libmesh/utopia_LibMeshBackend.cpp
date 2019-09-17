#include "utopia_LibMeshBackend.hpp"
#include "libmesh/petsc_vector.h"
#include "utopia_Adaptivity.hpp"
#include "libmesh/remote_elem.h"

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

       std::vector<SizeType> index;
       std::vector<SizeType> index_local;

       auto on_boundary = libMesh::MeshTools::find_boundary_nodes(V.mesh());

        std::cout << "apply_boundary_conditions Adaptivity begin: "  << std::endl;

      
      //  if(V.mesh().mesh_dimension()<3)
      //  {
      //       libMesh::MeshBase::const_node_iterator it = V.mesh().local_nodes_begin();
      //       const libMesh::MeshBase::const_node_iterator end_it = V.mesh().local_nodes_end();
      //       for ( ; it != end_it; ++it)
      //       {
      //           const libMesh::Node * node = *it;
                
      //           for (unsigned int comp = 0;comp < node->n_comp(V.equation_system().number(), 0); comp++)
      //           {
      //               const libMesh::dof_id_type node_dof = node->dof_number(V.equation_system().number(), 0, comp);
                    
      //               if(on_boundary.count(node->id()) && V.dof_map().is_constrained_dof(node_dof)) {

      //                    index.push_back(node_dof);
      //               }
      //           }
      //       }
      //   }

      // else

       {
            libMesh::MeshBase::const_element_iterator it = V.mesh().active_elements_begin();
            
            const libMesh::MeshBase::const_element_iterator end_it = V.mesh().active_elements_end();
            
            for ( ; it != end_it; ++it)
            {
                const libMesh::Elem * ele = *it;

                for(int kk=0; kk<ele->n_sides(); kk++)

                {             
                    auto side = ele->build_side_ptr(kk);

                    
                    if (ele->neighbor_ptr(kk) != libMesh::remote_elem) // V.mesh().boundary_info->boundary_ids(ele,kk).size
                    {

                        //std::cout<<"ciao, this is b_id=>"<<V.mesh().boundary_info->boundary_ids(ele,kk).at(0)<<std::endl;

                        index_local.clear();

                         //std::cout<<"side->n_nodes()"<<side->n_nodes()<<std::endl;

                        for (int ll=0; ll<ele->n_nodes(); ll++)
                        {

                           const libMesh::Node * node = ele->node_ptr(ll);


                           const libMesh::dof_id_type node_dof = node->dof_number(V.equation_system().number(), 0, 0);                

                            if(on_boundary.count(node->id()) && V.dof_map().is_constrained_dof(node_dof)) 
                            {
                                   
                                        index_local.push_back(node_dof);
                                        
                                        utopia::disp(node_dof);
           
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
                    

       

      


        const bool has_constaints = V.dof_map().constraint_rows_begin() != V.dof_map().constraint_rows_end();


        Size ls = local_size(mat);

        Size s = size(mat);

        set_zero_rows(mat, index, 1.);

        Write<UVector> w_v(vec, utopia::GLOBAL_INSERT);

        std::vector<SizeType> I(1,0);
        std::vector<double> value(1, 0);

        if(has_constaints) {
            libMesh::DofConstraintValueMap &rhs_values = V.dof_map().get_primal_constraint_values();

            Range r = range(vec);

            for(auto it=index.begin(); it < index.end(); ++it){
                int i = *it;
                std::cout<<"I=>"<<i<<std::endl;
                auto valpos = rhs_values.find(i);
                I[0] = i;
                value[0]=valpos->second;

                vec.set(I, value);

            }
        }
    

        c.stop();

        if(utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions end: " << c << std::endl;
        }
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


//for(auto it=index_local.begin(); it < index_local.end(); ++it)
//  {

//    //   auto check = (std::find(index.begin(), index.end(), *it) != index.end());

//    // if(!check)
//    // {
//      std::cout<<"ciao, this is node=>"<<*it<<std::endl;

//     //auto check = (std::find(index.begin(), index.end(), node_dof) != index.end());
//      index.push_back(*it);
//    // }
// }

