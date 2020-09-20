

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
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_ProjectedGaussSeidelNew.hpp"
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

namespace utopia {

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_tension_isotropic_phase_field_2(Input &in) {
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

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);
        PFFracTension2D<FunctionSpace> BC_setup(space);
        IncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace> > time_stepper(
            space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_tension_isotropic_phase_field_2);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_pressure_Tbar_isotropic_phase_field_2(Input &in) {
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

        InitialCondidtionPFTbar<FunctionSpace> IC_setup(space, 0.0);
        PFFracFixAllDisp<FunctionSpace> BC_setup(space);
        IncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace> > time_stepper(
            space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_pressure_Tbar_isotropic_phase_field_2);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_pressure_net_2_rmtr(Input &in) {
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
                             IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                             PFFracFixAllDisp<FunctionSpace>,
                             InitialCondidtionPFFracNet2D<FunctionSpace> >
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_pressure_net_2_rmtr);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_tension_phase_field_2_rmtr(Input &in) {
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
                             IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                             PFFracTension2D<FunctionSpace>,
                             InitialCondidtionPFTension<FunctionSpace> >
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_tension_phase_field_2_rmtr);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void AsphaltTension2d_rmtr(Input &in) {
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
                             IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                             AsphaltTension2D<FunctionSpace>,
                             AsphaltTension<FunctionSpace> >
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(AsphaltTension2d_rmtr);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void mixed2d_rmtr(Input &in) {
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
                             IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                             PFMixed2D<FunctionSpace>,
                             Mixed<FunctionSpace> >
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(mixed2d_rmtr);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void frac_plate_rmtr(Input &in) {
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
                             IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                             FracPlateBC<FunctionSpace>,
                             FracPlateIC<FunctionSpace> >
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(frac_plate_rmtr);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_pressure_network_isotropic_phase_field_2(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        // using Vector = typename FunctionSpace::Vector;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);

        stats.stop_and_collect("space-creation");
        stats.start();

        InitialCondidtionPFFracNet2D<FunctionSpace> IC_setup(space, 0.0);
        PFFracFixAllDisp<FunctionSpace> BC_setup(space);

        IncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace> > time_stepper(
            space, IC_setup, BC_setup);
        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_pressure_network_isotropic_phase_field_2);

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_isotropic_phase_field_3d(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);

        MLIncrementalLoading<FunctionSpace,
                             IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                             PFFracFixAllDisp<FunctionSpace>,
                             InitialCondidtionPFFracNet3D<FunctionSpace> >
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");
        space.comm().root_print(std::to_string(space.n_dofs()) + " coarse dofs");
        stats.stop_and_collect("space-creation");
    }

    UTOPIA_REGISTER_APP(petsc_isotropic_phase_field_3d);
}  // namespace utopia

#endif
