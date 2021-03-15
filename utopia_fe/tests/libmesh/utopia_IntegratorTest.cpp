#include "utopia_IntegratorTest.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_libmesh.hpp"

#include "utopia_BilinearFormExprIntegrator.hpp"
#include "utopia_CompositeEquationIntegrator.hpp"
#include "utopia_LinearFormExprIntegrator.hpp"

namespace utopia {

    void IntegratorTest::run(Input &in) {
        using ProductSpaceT = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;

        std::cout << "IntegratorTest" << std::endl;

        UIMesh<libMesh::DistributedMesh> mesh(this->comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));

        in.get("mesh", mesh);
        in.get("space", space);

        auto &V = space.space();

        std::shared_ptr<BilinearIntegrator<ProductSpaceT>> bilinear_integrator;
        std::shared_ptr<EquationIntegrator<ProductSpaceT>> eq_integrator;
        std::shared_ptr<LinearIntegrator<ProductSpaceT>> linear_integrator;

        linear_integrator = linear_form(
            space.space_ptr(), surface_integral(inner(coeff(1.0), test(V[0])) + inner(coeff(1.0), test(V[1]))));

        UVector vec;
        utopia::assemble(*linear_integrator, vec);

        double vol = sum(vec);
        approxeq(8.0, vol, 1e-7);

        bilinear_integrator = bilinear_form(
            space.space_ptr(), surface_integral(inner(trial(V[0]), test(V[0])) + inner(trial(V[1]), test(V[1]))));

        USparseMatrix mat;
        utopia::assemble(*bilinear_integrator, mat);
        double vol_mat = sum(mat);

        approxeq(8.0, vol_mat, 1e-7);

        eq_integrator = equation(bilinear_integrator, linear_integrator);

        utopia::assemble(*eq_integrator, mat, vec);
        vol = sum(vec);

        approxeq(8.0, vol, 1e-7);

        vol_mat = sum(mat);

        approxeq(8.0, vol_mat, 1e-7);

        CompositeEquationIntegrator<ProductSpaceT> composite_integrator;

        composite_integrator.add_integrator(eq_integrator);
        composite_integrator.add_integrator(eq_integrator);

        utopia::assemble(composite_integrator, mat, vec);

        vol = sum(vec);

        approxeq(16.0, vol, 1e-7);

        vol_mat = sum(mat);

        approxeq(16.0, vol_mat, 1e-7);
    }
}  // namespace utopia
