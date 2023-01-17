

#include "utopia_Base.hpp"

// FIXME: this file causes nvcc to fail
#ifndef KOKKOS_ENABLE_CUDA

#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_BCSetup.hpp"
#include "utopia_Backtracking.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_FracNetGenerator2D.hpp"
#include "utopia_FracNetGenerator3D.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_IncrementalLoading.hpp"
#include "utopia_InitialCondition.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_LBFGS.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MLIncrementalLoading.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PhaseFieldVolDevSplit.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"

#include "utopia_QuasiNewtonBound.hpp"
#include "utopia_QuasiTrustRegionVariableBound.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_TrustRegionVariableBound.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_RedundantQPSolver.hpp"
#endif  // UTOPIA_WITH_PETSC

#include <chrono>
#include <cmath>
#include <random>

#ifdef UTOPIA_WITH_VC
#include "utopia_vc_IsotropicPhaseFieldForBrittleFractures.hpp"
#endif  // UTOPIA_WITH_VC

namespace utopia {

    template <class FunctionSpace>
#ifdef UTOPIA_WITH_VC
    using FractureModel = utopia::VcIsotropicPhaseFieldForBrittleFractures<FunctionSpace>;
#else
    using FractureModel = utopia::IsotropicPhaseFieldForBrittleFractures<FunctionSpace>;
#endif

    // using FractureModel = utopia::PhaseFieldVolDevSplit<FunctionSpace>;

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void franetg(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        stats.start();

        MLIncrementalLoading<FunctionSpace,
                             FractureModel<FunctionSpace>,
                             UniaxialLoading2D<FunctionSpace>,
                             UniformSpacingOnLine<FunctionSpace>>
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(franetg);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void VolDevSplit(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        stats.start();

        MLIncrementalLoading<FunctionSpace,
                             utopia::PhaseFieldVolDevSplit<FunctionSpace>,
                             BiaxialLoading2D<FunctionSpace>,
                             SlantedCrack2D<FunctionSpace>>
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(VolDevSplit);


}  // namespace utopia

#endif
