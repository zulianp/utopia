


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

#include "utopia_LinearElasticityFE.hpp"

#include "utopia_app_utils.hpp"

#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"

#include <cmath>

namespace utopia {

    template<class FunctionSpace>
    void neumman_test(FunctionSpace &space, Input &in)
    {
        using Point  = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        NeumannBoundaryCondition<FunctionSpace> bc(
            space,
            //selector
            SideSet::left(),
            //value
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 1.0;
            },
            0
        );

        Vector v;
        space.create_vector(v);

        bc.apply(v);

        Scalar side_area = sum(v);

        space.comm().root_print(side_area);

        space.write("neumman.vtr", v);

        NeumannBoundaryCondition<FunctionSpace> bc2(
            space,
            //selector
            UTOPIA_LAMBDA(const Point &x) -> bool {
                return x[0] <= device::epsilon<Scalar>() && (x[1] >= 0.5 && x[1] <= 1.0);
            },
            //value
            UTOPIA_LAMBDA(const Point &x) -> Scalar {
                return x[1];
            },
            0
        );

        v.set(0.0);
        bc2.apply(v);

        side_area = sum(v);

        space.comm().root_print(side_area);

        bc2.apply(v);

        side_area = sum(v);

        space.comm().root_print(side_area);

        space.write("neumman2.vtr", v);
    }

    void neumann_example_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = 1;

        using Mesh          = utopia::PetscDM<Dim>;
        using Elem          = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;

        FunctionSpace space;
        space.read(in);

        neumman_test(space, in);
    }

    UTOPIA_REGISTER_APP(neumann_example_2);

    void neumann_example_3(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = 1;

        using Mesh          = utopia::PetscDM<Dim>;
        using Elem          = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;

        FunctionSpace space;
        space.read(in);

        neumman_test(space, in);
    }


    UTOPIA_REGISTER_APP(neumann_example_3);
}
