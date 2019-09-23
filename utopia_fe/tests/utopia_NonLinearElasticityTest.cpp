#include "utopia_libmesh.hpp"
#include "utopia_NonLinearElasticityTest.hpp"

#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include <algorithm>

namespace utopia {

    std::shared_ptr<libMesh::DistributedMesh> square(
        libMesh::Parallel::Communicator &comm,
        const unsigned int nx,
        const unsigned int ny,
        const double x_left,
        const double x_right,
        const double y_bottom,
        const double y_top,
        const libMesh::ElemType &elem_type)
    {
        auto mesh = std::make_shared<libMesh::DistributedMesh>(comm);

        libMesh::MeshTools::Generation::build_square(*mesh,
            nx, ny,
            x_left, x_right,
            y_bottom, y_top,
            elem_type
            );

        return mesh;
    }

        //saint-venant kirchoff
        // auto C = F_t * F;
        // auto E = 0.5 * (C - identity());
        // auto S = 2.0 * mu * E + lambda * (trace(E) * identity());
        // auto P = F * S;

        // auto strain_lin = 0.5 * (F_t * grad(u) + transpose(grad(u)) * F);
        // auto stress_lin = F * (
        // 	2.0 * mu * strain_lin + lambda * (trace(strain_lin) * identity())
        // 	) + grad(u) * S;

        /////////////////////////////




    void NonLinearElasticityTest::run(Input &in)
    {
        assert(false);
        // auto mesh = square(this->comm(), 20, 20, 0., 1., 0., 1., libMesh::QUAD8);

        // const auto dim = mesh->mesh_dimension();

        // auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
        // auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("neo-hookean");

        // const double mu = 1.;
        // const double lambda = 1.;

        // auto elem_order = libMesh::SECOND;

        // ////////////////////////////////////////////

        // auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
        // auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
        // auto V = Vx * Vy;

        // if(dim == 3) {
        //     V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
        // }

        // auto u = trial(V);
        // auto v = test(V);

        // auto ux = u[0];
        // auto uy = u[1];

        // ////////////////////////////////////////////
        // UVector sol;
        // auto uk = interpolate(sol, u);

        // auto F 		 = identity() + grad(uk);
        // auto F_t 	 = transpose(F);
        // auto F_inv   = inv(F);
        // auto F_inv_t = transpose(F_inv);
        // auto J       = det(F);

        // //compressible neo-hookean
        // auto P = mu * (F - F_inv_t) + (lambda * logn(J)) * F_inv_t;

        // auto stress_lin = mu * grad(u)
        // -(lambda * logn(J) - mu) * F_inv_t * transpose(grad(u)) * F_inv_t
        // + inner(lambda * F_inv_t, grad(u)) * F_inv_t;

        // auto l_form = inner(P, grad(v)) * dX;
        // auto b_form = inner(stress_lin, grad(v)) * dX;
        // auto eq = b_form == l_form;

        // ////////////////////////////////////////////

        // libMesh::Nemesis_IO io(*mesh);

        // if(dim == 3) {
        //     auto uz = u[2];

        //     nl_solve(
        //         equations(
        //                 eq
        //             ),
        //         constraints(
        //             boundary_conditions(ux == coeff(0.),    {1, 3}),
        //             boundary_conditions(uz == coeff(0.),    {1, 3}),
        //             boundary_conditions(uy == coeff(0.05),  {1}),
        //             boundary_conditions(uy == coeff(-0.05), {3})
        //             ),
        //         sol
        //         );
        // } else {

        //     nl_solve(
        //         equations(
        //              eq
        //             ),
        //         constraints(
        //             boundary_conditions(ux == coeff(0.),  {0, 2}),
        //             boundary_conditions(uy == coeff(0.1), {0}),
        //             boundary_conditions(uy == coeff(-0.1), {2})
        //             ),
        //         sol
        //         );
        // }



        // convert(sol, *sys.solution);
        // sys.solution->close();
        // io.write_equation_systems("hyper-elasticity.e", *equation_systems);
    }
}