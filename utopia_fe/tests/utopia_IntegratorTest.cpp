#include "utopia_IntegratorTest.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIFunctionSpace.hpp"

#include "utopia_LinearFormExprIntegrator.hpp"
#include "utopia_BilinearFormExprIntegrator.hpp"

namespace utopia {

    void IntegratorTest::run(Input &in)
    {
        using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;

        std::cout << "IntegratorTest" << std::endl;

        UIMesh<libMesh::DistributedMesh> mesh(this->comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));

        in.get("mesh", mesh);
        in.get("space", space);

        auto &V = space.space();

        auto linear_integrator = linear_form(
            space.space_ptr(),
            surface_integral(
                        inner(coeff(1.0), test(V[0])) + 
                        inner(coeff(1.0), test(V[1]))
                    )
        );

        UVector vec;
        utopia::assemble(*linear_integrator, vec);

        double vol = sum(vec);
        disp(vol);

        auto bilinear_integrator = bilinear_form(
            space.space_ptr(), 
            surface_integral(
                        inner(trial(V[0]), test(V[0])) + 
                        inner(trial(V[1]), test(V[1]))
                    )
        );

        USparseMatrix mat;
        utopia::assemble(*bilinear_integrator, mat);
        double vol_mat = sum(mat);

        disp(vol_mat);
    }
}

