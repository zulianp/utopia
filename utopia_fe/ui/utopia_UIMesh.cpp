#include "utopia_UIMesh.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_TransferAssembler.hpp"

#include "moonolith_communicator.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"


namespace utopia {
    
    void refine_at_intersection(
        const std::shared_ptr<libMesh::UnstructuredMesh> &fracture_network,
        const libMesh::Order &elem_order,
        const std::shared_ptr<libMesh::UnstructuredMesh> &mesh,
        const int refinement_loops,
        const bool use_interpolation
        )
    {


        
        libMesh::MeshRefinement mesh_refinement(*mesh);


        for(int i = 0; i < refinement_loops; ++i) {
        //equations system
            auto vol_equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
            auto &vol_sys = vol_equation_systems->add_system<libMesh::LinearImplicitSystem>("vol_sys");

            auto surf_equation_systems = std::make_shared<libMesh::EquationSystems>(*fracture_network);
            auto &surf_sys = surf_equation_systems->add_system<libMesh::LinearImplicitSystem>("surf_sys");

        //scalar function space
            auto V_vol  = LibMeshFunctionSpace(vol_equation_systems, libMesh::LAGRANGE, elem_order,      "u_vol");
            auto V_surf = LibMeshFunctionSpace(surf_equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_surf");

            V_vol.initialize();
            V_surf.initialize();

            Chrono c;
            c.start();
            USparseMatrix B;

            // mesh->print_info();
            // fracture_network->print_info();
            moonolith::Communicator comm(mesh->comm().get());
            if(assemble_volume_transfer(
                comm,
                mesh,
                fracture_network,
                make_ref(V_vol.dof_map()),
                make_ref(V_surf.dof_map()),
                0,
                0,
                true,
                1,
                B,
                {},
                use_interpolation))
            {
                c.stop();
                std::cout << c << std::endl;

                Interpolator interp(make_ref(B));
                interp.normalize_rows();
                interp.describe(std::cout);

                USparseMatrix D_inv = diag(1./sum(B, 1));
                USparseMatrix T = D_inv * B;

                USparseMatrix T_t = transpose(T);
                UVector t_temp = sum(T_t, 1);
                UVector t = ghosted(local_size(t_temp).get(0), size(t_temp).get(0), V_vol.dof_map().get_send_list());
                t = t_temp;


                std::vector<libMesh::dof_id_type> indices;
                std::vector<double> values;

                mesh_refinement.clean_refinement_flags();

                Read<UVector> r_(t);
                for(auto e_it = elements_begin(*mesh); e_it != elements_end(*mesh); ++e_it) {
                    V_vol.dof_map().dof_indices(*e_it, indices);
                    t.get(indices, values);

                    double val = std::accumulate(values.begin(), values.end(), 0.,  std::plus<double>());
                    if(val > 0) {
                        (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
                    }

                }

                mesh_refinement.make_flags_parallel_consistent();
                mesh_refinement.refine_elements();
                mesh_refinement.test_level_one(true);

            } else {
                assert(false);
            }
        }

        // mesh_refinement.clean_refinement_flags();
        mesh->prepare_for_use();
    }

}

