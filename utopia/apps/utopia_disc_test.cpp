#include "utopia_Base.hpp"

//include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Core.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_petsc.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_PhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"

#include <cmath>

namespace utopia {

    void disc_test(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = 1;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using Point            = Mesh::Point;
        using Scalar           = Mesh::Scalar;

        FunctionSpace space;
        space.read(in);

        //f(w) = 4*x^3 + y^2 + 3
        //lapl f(w) = 4*6*x + 2

        for(auto s : space.mesh().sides()) {
            space.emplace_dirichlet_condition(
                s,
                UTOPIA_LAMBDA(const Point &x) -> Scalar {
                    return 4 * x[0] * x[0] * x[0] + x[1] * x[1] + 3;
                },
                0
            );
        }

        PoissonFE<FunctionSpace> problem(space);

        problem.init_forcing_function(UTOPIA_LAMBDA(const Point &x) -> Scalar {
            return (4.0 * 6.0) * x[0] + 2.0;
        });


        PetscVector x, g;
        PetscMatrix H;

        space.create_vector(x);

        problem.update(x);

        problem.gradient(x, g);
        problem.hessian(x, H);

        space.apply_constraints(H, g);

        solve(H, g, x);
        space.write("X.vts", x);
    }

    UTOPIA_REGISTER_APP(disc_test);

}