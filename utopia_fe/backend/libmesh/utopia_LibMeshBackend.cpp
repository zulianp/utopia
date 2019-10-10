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

       std::vector<int> index, dirichel_id;
       std::vector<SizeType> index_local;


       Adaptivity::compute_boundary_nodes(V.mesh(), V.dof_map(), V.equation_system().number(), 0, index);
  

        const bool has_constaints = V.dof_map().constraint_rows_begin() != V.dof_map().constraint_rows_end();


        Size ls = local_size(mat);

        Size s = size(mat);

        set_zero_rows(mat, index, 1.);

        Write<UVector> w_v(vec, utopia::GLOBAL_INSERT);

        std::vector<SizeType> I(1,0);
        std::vector<double> value(1, 0);

        if(has_constaints) 
        {
            libMesh::DofConstraintValueMap &rhs_values = V.dof_map().get_primal_constraint_values();

            Range r = range(vec);

            for(auto it=index.begin(); it < index.end(); ++it)
            {
                int i = *it;
                auto valpos = rhs_values.find(i);
                I[0] = i;
                value[0]=valpos->second;


                // if (V.mesh().processor_id()==0) std::cout<<"i"<<i<<"=>"<<value[0]<<std::endl;
                // if (V.mesh().processor_id()==1) std::cout<<"i"<<i<<"=>"<<value[0]<<std::endl;
                // if (V.mesh().processor_id()==2) std::cout<<"i"<<i<<"=>"<<value[0]<<std::endl;
                // if (V.mesh().processor_id()==3) std::cout<<"i"<<i<<"=>"<<value[0]<<std::endl;
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



