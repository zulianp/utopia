
#include "utopia_fe_kokkos_fix.hpp"

#include "utopia_BoundaryIntegralTest.hpp"

#include "utopia.hpp"

// fe extension
#include "moonolith_profiler.hpp"
#include "utopia_ContactSimParams.hpp"
// #include "utopia_Socket.hpp"
#include "utopia_fe_EDSL.hpp"

#include <libmesh/const_function.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/parallel_mesh.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/petsc_vector.h>
#include "libmesh/linear_partitioner.h"

#include "utopia_libmesh_NonLinearFEFunction.hpp"

#include "libmesh/mesh_base.h"
#include "libmesh/nemesis_io.h"

#include <iostream>

namespace utopia {

    void BoundaryIntegralTest::run(Input &in) {
        moonolith::Communicator comm(this->comm().get());

        const unsigned int nx = 6;
        const unsigned int ny = 6;
        const unsigned int nz = 6;

        auto elem_order = libMesh::SECOND;
        auto mesh = std::make_shared<libMesh::DistributedMesh>(this->comm());

        libMesh::MeshTools::Generation::build_cube(*mesh, nx, ny, nz, 0, 1, 0, 1., 0, 1., libMesh::TET4);

        mesh->all_second_order(true);

        auto dim = mesh->mesh_dimension();

        auto sys = std::make_shared<libMesh::EquationSystems>(*mesh);
        sys->add_system<libMesh::LinearImplicitSystem>("bit");

        auto V = LibMeshFunctionSpace(sys, libMesh::LAGRANGE, elem_order, "u");
        V.initialize();

        auto u = trial(V);
        auto v = test(V);

        libMesh::QGauss q_gauss(dim - 1, libMesh::FOURTH);
        q_gauss.init(libMesh::TRI6);

        auto fe = libMesh::FEBase::build(mesh->mesh_dimension(), elem_order);
        fe->get_phi();
        fe->get_JxW();

        auto &dof_map = V.dof_map();

        std::vector<libMesh::dof_id_type> dof_indices;

        auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
                                  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

        USparseMatrix boundary_mass_matrix = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

        double mat_sum = 0.0;

        {
            Write<USparseMatrix> w_b(boundary_mass_matrix);
            LMDenseMatrix elemat;

            // for(const auto & elem : mesh->active_local_element_ptr_range()) {
            for (auto e_it = mesh->active_local_elements_begin(); e_it != mesh->active_local_elements_end(); ++e_it) {
                auto elem = *e_it;

                dof_map.dof_indices(elem, dof_indices);
                elemat.resize(dof_indices.size(), dof_indices.size());
                elemat.set(0.0);

                bool has_assembled = false;

                // for(auto side : elem->side_index_range()) {
                for (auto side = 0; side < elem->n_sides(); ++side) {
                    if ((elem->neighbor_ptr(side) != libmesh_nullptr)) {
                        continue;
                    }

                    fe->attach_quadrature_rule(&q_gauss);
                    fe->reinit(elem, side);

                    auto &phi = fe->get_phi();
                    auto &JxW = fe->get_JxW();

                    for (std::size_t i = 0; i < dof_indices.size(); ++i) {
                        for (std::size_t j = 0; j < dof_indices.size(); ++j) {
                            for (std::size_t qp = 0; qp < JxW.size(); ++qp) {
                                elemat.add(i, j, phi[i][qp] * phi[j][qp] * JxW[qp]);
                            }
                        }
                    }

                    has_assembled = true;
                }

                if (has_assembled) {
                    add_matrix(elemat, dof_indices, dof_indices, boundary_mass_matrix);
                    mat_sum += elemat.sum();
                }
            }
        }

        USparseMatrix boundary_mass_matrix_2;
        assemble(inner(u, v) * dS, boundary_mass_matrix_2);

        double surface_area = sum(boundary_mass_matrix);
        double surface_area_2 = sum(boundary_mass_matrix_2);

        utopia_test_assert(approxeq(surface_area, 6., 1e-10));
        utopia_test_assert(approxeq(surface_area_2, 6., 1e-10));

        USparseMatrix diff_mm = boundary_mass_matrix - boundary_mass_matrix_2;
        const double diff_norm = norm2(diff_mm);

        utopia_test_assert(approxeq(diff_norm, 0.));

        USparseMatrix side_mass_matrix;
        assemble(surface_integral(inner(u, v), 1), side_mass_matrix);

        UVector side_mass_vec = sum(side_mass_matrix, 1);

        const double side_area_1 = sum(side_mass_matrix);
        utopia_test_assert(approxeq(side_area_1, 1., 1e-10));

        USparseMatrix side_mass_matrix_12;
        assemble(surface_integral(inner(u, v), 1) + surface_integral(inner(u, v), 2), side_mass_matrix_12);
        const double side_area_12 = sum(side_mass_matrix_12);
        utopia_test_assert(approxeq(side_area_12, 2., 1e-10));

        USparseMatrix boundary_lapl;
        assemble(surface_integral(inner(grad(u), grad(v))), boundary_lapl);
        const double sum_lapl = sum(boundary_lapl);
        utopia_test_assert(approxeq(sum_lapl, 0., 1e-10));

        UVector b_fun;
        assemble(surface_integral(inner(coeff(1.), v), 1), b_fun);
    }
}  // namespace utopia
