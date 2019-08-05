#include "utopia_libmesh.hpp"
#include "utopia_MechTest.hpp"

#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"


#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_Mechanics.hpp"
#include "utopia_AffineTransform.hpp"
#include "utopia_Contact.hpp"
#include "utopia_LinearElasticity.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

#include <algorithm>
#include <memory>
#include <array>

typedef std::array<double, 2> Point2d;
typedef std::array<double, 3> Point3d;
typedef std::function<std::array<double, 2>(const std::array<double, 2> &p)> Fun2d;

namespace utopia {

    void MechTest::run(Input &in)
    {
        auto mesh = std::make_shared<libMesh::DistributedMesh>(this->comm());

        // const unsigned int n = 8;
        const unsigned int dim = 2;
        const double dt = 1;

        mesh->read("../data/ironing.e");

        // libMesh::MeshTools::Generation::build_square(
        // 	*mesh,
        // 	n, n,
        // 	0, 1,
        // 	0, 1.,
        // 	libMesh::QUAD8);

        auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
        auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("mech_test");

        auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST,  "disp_x");
        auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "disp_y");
        auto V  = Vx * Vy;

        auto u  = trial(V);
        auto ux = trial(Vx);
        auto uy = trial(Vy);

        double t = 0.05;
        double angle = (60./180*M_PI);
        Fun2d trafo = [&t, angle](const Point2d &p) -> Point2d {
            AffineTransform trafo;
            trafo.make_rotation(2, t*angle, 'y');
            trafo.translation[0] = -p[0];
            trafo.translation[1] = -p[1];
            return trafo.apply(p);
        };

        auto constr = constraints(
            boundary_conditions(ux == coeff(0.), {3}),
            boundary_conditions(uy == coeff(-0.05), {3}),
            boundary_conditions(u  == coeff(trafo), {4})
        );

        FEBackend<LIBMESH_TAG>::init_constraints(constr);
        Vx.initialize();

        Size gs({Vx.dof_map().n_dofs()});
        Size ls({Vx.dof_map().n_local_dofs()});

        UVector internal_force;
        USparseMatrix stiffness_matrix;

        MechanicsContext mech_ctx;
        MechanicsState old;
        MechanicsState current;

        mech_ctx.init_mass_matrix(V);

        old.init(ls, gs);
        current.init(ls, gs);

        LameeParameters params;
        auto elast = std::make_shared<LinearElasticity<decltype(V), USparseMatrix, UVector>>(V, params);
        elast->assemble_hessian_and_gradient(old.displacement, mech_ctx.stiffness_matrix, old.internal_force);

        auto vx = test(Vx);
        auto vy = test(Vy);

        auto ef = std::make_shared<ConstantExternalForce>();
        // ef->init((inner(coeff(0.5), vx) + inner(coeff(-0.1), vy)) * dX);
        ef->init((inner(coeff(0.), vx) + inner(coeff(0.), vy)) * dX);

        apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, ef->value);
        ef->eval(old.t, old.external_force);
        ef->eval(old.t + dt, current.external_force);

        auto integrator = std::make_shared<ImplicitEuler>(dim, Vx.dof_map());

        Contact contact;
        ContactParams contact_params;
        contact_params.contact_pair_tags = {{2, 1}};
        contact_params.search_radius = 0.2;
        if(!contact.init(mesh, make_ref(Vx.dof_map()), contact_params)) {
            std::cerr << "[Error] contact failed" << std::endl;
        }

        integrator->apply(dt, mech_ctx, contact, Friction(), old, current);

        convert(current.displacement, *sys.solution);
        sys.solution->close();
        libMesh::ExodusII_IO(*mesh).write_equation_systems ("mech_test.e", *equation_systems);
    }
}
