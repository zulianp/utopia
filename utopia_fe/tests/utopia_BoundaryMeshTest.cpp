#include "utopia_BoundaryMeshTest.hpp"
#include "libmesh/boundary_mesh.h"
#include "moonolith_synched_describable.hpp"
#include "utopia_libmesh.hpp"

namespace utopia {
    class BoundaryMeshTestImpl {
    public:
        using FunctionSpaceT = utopia::LibMeshFunctionSpace;

        BoundaryMeshTestImpl(libMesh::Parallel::Communicator &comm) : comm(comm), n(2) { build(); }

        void run() {
            // UTOPIA_UNIT_TEST_BEGIN("BounaryMeshTest");
            UTOPIA_RUN_TEST(dof_relation);
            UTOPIA_RUN_TEST(mass_matrix);
            // UTOPIA_UNIT_TEST_END("BounaryMeshTest");
        }

        void dof_relation() {
            auto V = FunctionSpaceT(*mesh, libMesh::LAGRANGE, libMesh::FIRST, "u");
            V.initialize();

            libMesh::BoundaryMesh b_mesh(comm, mesh->mesh_dimension() - 1);
            mesh->boundary_info->sync(b_mesh);

            auto V_surf = FunctionSpaceT(b_mesh, libMesh::LAGRANGE, libMesh::FIRST, "u");
            V_surf.initialize();

            auto &dof_map = V.dof_map();
            auto &dof_map_surf = V_surf.dof_map();

            std::vector<libMesh::dof_id_type> side_indices, indices, surf_indices;

            std::stringstream ss;
            for (auto e_it = elements_begin(*mesh); e_it != elements_end(*mesh); ++e_it) {
                auto &e = **e_it;

                auto n_sides = e.n_sides();

                dof_map.dof_indices(&e, indices);

                // ss << "-----------------------------\n";
                // ss << "-----------------------------\n";
                // print_vector(std::begin(indices), std::end(indices), ss);
                // ss << "....................................\n";

                for (std::size_t i = 0; i < n_sides; ++i) {
                    auto side = e.build_side_ptr(i);

                    dof_map.dof_indices(side.get(), side_indices);

                    // print_vector(std::begin(side_indices), std::end(side_indices), ss);

                    dof_map_surf.dof_indices(side.get(), surf_indices);

                    // print_vector(std::begin(surf_indices), std::end(surf_indices), ss);
                    // ss << "....................................\n";

                    utopia_test_assert(side_indices == surf_indices);
                }

                // ss << "-----------------------------\n";
            }

            moonolith::Communicator m_comm(comm.get());
            // moonolith::synch_describe(ss.str(), m_comm, std::cout);
        }

        void mass_matrix() {
            auto V_vol = FunctionSpaceT(*mesh, libMesh::LAGRANGE, libMesh::FIRST, "u");
            V_vol.initialize();

            auto form_vol = inner(trial(V_vol), test(V_vol)) * dX;

            USparseMatrix mat_vol;
            utopia::assemble(form_vol, mat_vol);

            double volume = sum(mat_vol);
            utopia_test_assert(approxeq(volume, 1.));

            libMesh::BoundaryMesh b_mesh(comm, mesh->mesh_dimension() - 1);
            mesh->boundary_info->sync(b_mesh);

            auto V = FunctionSpaceT(b_mesh, libMesh::LAGRANGE, libMesh::FIRST, "u");
            V.initialize();

            auto form = inner(trial(V), test(V)) * dX;

            USparseMatrix mat;
            utopia::assemble(form, mat);

            double surface_area = sum(mat);
            utopia_test_assert(approxeq(surface_area, (mesh->mesh_dimension() == 2 ? 4. : 6.)));
        }

        void build() {
            mesh = std::make_shared<libMesh::DistributedMesh>(comm);
            // libMesh::MeshTools::Generation::build_square(
            // 	*mesh,
            // 	n, n,
            // 	0, 1,
            // 	0, 1.,
            // 	libMesh::QUAD4
            // );

            libMesh::MeshTools::Generation::build_cube(*mesh, n, n, n, 0, 1, 0, 1., 0, 1., libMesh::HEX8);
        }

    private:
        libMesh::Parallel::Communicator &comm;
        std::shared_ptr<libMesh::UnstructuredMesh> mesh;
        int n;
    };

    void BoundaryMeshTest::run(Input &in) {
        BoundaryMeshTestImpl impl(comm());
        impl.run();
    }
}  // namespace utopia
