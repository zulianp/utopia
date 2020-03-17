


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

#include <cmath>

namespace utopia {

    void neumann_example(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using Point            = FunctionSpace::Point;
        using Scalar           = FunctionSpace::Scalar;

        FunctionSpace space;
        space.read(in);

        NeumannBoundaryCondition<FunctionSpace> bc(space);
    }

    UTOPIA_REGISTER_APP(neumann_example);
}
