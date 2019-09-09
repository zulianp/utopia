#include "utopia_LibMeshBackend.hpp"
#include "libmesh/petsc_vector.h"
#include "utopia_Adaptivity.hpp"

namespace utopia {

    void apply_boundary_conditions(LibMeshFunctionSpace &V, USparseMatrix &mat, UVector &vec)
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

       // V.mesh().get_boundary_info().print_info();

       // std::vector<libMesh::boundary_id_type> node_boundaries;

       // V.mesh().get_boundary_info().build_node_boundary_ids(node_boundaries);

       // std::vector<libMesh::dof_id_type> node_id_list;
       

       auto on_boundary = libMesh::MeshTools::find_boundary_nodes(V.mesh());


       {
            libMesh::MeshBase::const_node_iterator it = V.mesh().local_nodes_begin();
            const libMesh::MeshBase::const_node_iterator end_it = V.mesh().local_nodes_end();
            for ( ; it != end_it; ++it)
            {
                const libMesh::Node * node = *it;
                
                for (unsigned int comp = 0;comp < node->n_comp(V.equation_system().number(), 0); comp++)
                {
                    const libMesh::dof_id_type node_dof = node->dof_number(V.equation_system().number(), 0, comp);

                     //std::cout<<"node_dof "<<node_dof <<std::endl;
                    
                    if(on_boundary.count(node->id()) && V.dof_map().is_constrained_dof(node_dof)) {
                        
                        index.push_back(node_dof);
                    }                   
                }
            }
       }


        const bool has_constaints = V.dof_map().constraint_rows_begin() != V.dof_map().constraint_rows_end();


        Size ls = local_size(mat);

        Size s = size(mat);

        set_zero_rows(mat, index, 1.);

        Write<UVector> w_v(vec);

        if(has_constaints) {
            libMesh::DofConstraintValueMap &rhs_values = V.dof_map().get_primal_constraint_values();

            Range r = range(vec);

            for(auto it=index.begin(); it < index.end(); ++it){
                int i = *it;
                auto valpos = rhs_values.find(i);
                vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);

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

