

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

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_GenericPhaseFieldFormulation.hpp"
#include "utopia_IsotropicGenericPhaseField.hpp"
#include "utopia_PhaseFieldDerivativeCheck_Original.hpp"
#include "utopia_TensionSplitGenericPhaseField.hpp"
#include "utopia_VolDevGenericPhaseField.hpp"
#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_RedundantQPSolver.hpp"
#endif  // UTOPIA_ENABLE_PETSC

#include <chrono>
#include <cmath>
#include <random>

#ifdef UTOPIA_ENABLE_VC
#include "utopia_vc_IsotropicPhaseFieldForBrittleFractures.hpp"
#endif  // UTOPIA_ENABLE_VC

namespace utopia {

    template <class FunctionSpace>
#ifdef UTOPIA_ENABLE_VC
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

        //        IsotropicGenericPhaseField<FunctionSpace, Dim, AT2> BrittleRock(space);

        MLIncrementalLoading<FunctionSpace,
                             // FractureModel<FunctionSpace>,
                             IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>,
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
    static void HomogeneousBarPseudo1D(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting HomogeneousBarPseudo1D Model" << std::endl;

        stats.start();

        MLIncrementalLoading<FunctionSpace,
                             utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>,
                             UniaxialLoading2D<FunctionSpace>,
                             HomogeneousBar<FunctionSpace>>
            time_stepper(space);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(HomogeneousBarPseudo1D);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void HomogeneousBarPseudo1DSingleLevel(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        //         using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>; //works
        //        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT2_CUBIC>; //works
        //        using ProblemType = utopia::IsotropicPhaseFieldForBrittleFractures<FunctionSpace, Dim>; //works
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;
        //        using ProblemType = utopia::PhaseFieldVolDevSplit<FunctionSpace, Dim>;
        //        using ProblemType = utopia::TensionSplitGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting HomogeneousBarPseudo1DSingleLevel Model" << std::endl;

        stats.start();

        HomogeneousBar<FunctionSpace> IC_setup(space, 0.0);
        UniaxialLoading2D<FunctionSpace> BC_setup(space, 0.0);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(HomogeneousBarPseudo1DSingleLevel);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void BoxForceTest(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        //         using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>; //works
        //        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT2_CUBIC>; //works
        //        using ProblemType = utopia::IsotropicPhaseFieldForBrittleFractures<FunctionSpace, Dim>; //works
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;
        //        using ProblemType = utopia::PhaseFieldVolDevSplit<FunctionSpace, Dim>;
        //        using ProblemType = utopia::TensionSplitGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting HomogeneousBarPseudo1DSingleLevel Model" << std::endl;

        stats.start();

        HomogeneousBar<FunctionSpace> IC_setup(space, 0.0);

        DirichletAndVolConstraints<FunctionSpace, BoxLoading<FunctionSpace>, NoDamage<FunctionSpace>> BC_setup(space);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(BoxForceTest);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void HomogeneousBarDerivativeCheck(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        //        using ProblemType = utopia::PhaseFieldDerivativeCheck_Original<FunctionSpace, Dim>;
        using ProblemType = utopia::IsotropicPhaseFieldForBrittleFractures<FunctionSpace, Dim>;
        //        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>; //Fixed
        //        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT2>; //Fixed
        //        using ProblemType = utopia::PhaseFieldVolDevSplit<FunctionSpace, Dim>; //Fixed
        //        using ProblemType = utopia::TensionSplitGenericPhaseField<FunctionSpace, Dim, AT1>; //Checking

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting HomogeneousBarDerivativeCheck Model" << std::endl;

        stats.start();

        //        HomogeneousBar<FunctionSpace> IC_setup(space, 0.0);
        SmoothQuadraticPhaseField<FunctionSpace> IC_setup(space, 0.0);
        BiaxialLoading2D<FunctionSpace> BC_setup(space, 0.0);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(HomogeneousBarDerivativeCheck);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void TensionTestSingleLevel(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>;
        //        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting TensionTest SL Model" << std::endl;

        stats.start();

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);
        PFFracTension<FunctionSpace> BC_setup(space, 0.0);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(TensionTestSingleLevel);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void TensionTestSingleLevel_3D(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting TensionTest SL Model" << std::endl;

        stats.start();

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);
        PFFracTension<FunctionSpace> BC_setup(space, 0.0);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(TensionTestSingleLevel_3D);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void ShearTestSingleLevel(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        //        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT2>;
        //        using ProblemType = utopia::IsotropicPhaseFieldForBrittleFractures<FunctionSpace, Dim>;
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting ShearTest SL Model ISO" << std::endl;

        stats.start();

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);
        PFFracShear2D<FunctionSpace> BC_setup(space, 0.0);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(ShearTestSingleLevel);

    // // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // // //
    static void FaultTest(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        //        using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT2>;
        //        using ProblemType = utopia::IsotropicPhaseFieldForBrittleFractures<FunctionSpace, Dim>;
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting ShearTest SL Model ISO" << std::endl;

        stats.start();

        SingleFault<FunctionSpace> IC_setup(space, 0.0);
        PFFracShear2D<FunctionSpace> BC_setup(space, 0.0);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(FaultTest);

    static void SingleSedimentaryLayer(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        // using ProblemType = utopia::TensionSplitGenericPhaseField<FunctionSpace, Dim, AT1>;
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;
        // using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting Single Layer Model" << std::endl;

        stats.start();

        DamagedSedimentaryLayers<FunctionSpace> IC_setup(space, 0.0);
        DirichletAndVolConstraints<FunctionSpace, SedimentaryLayers_BC<FunctionSpace>, LayeredSubdomain<FunctionSpace>>
            BC_setup(space);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(SingleSedimentaryLayer);

    static void SingleSedimentaryLayer_Reg(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;
        // using ProblemType = utopia::TensionSplitGenericPhaseField<FunctionSpace, Dim, AT1>;
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1_Regularised>;
        // using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting Single Layer Model" << std::endl;

        stats.start();

        DamagedSedimentaryLayers<FunctionSpace> IC_setup(space, 0.0);
        DirichletAndVolConstraints<FunctionSpace, SedimentaryLayers_BC<FunctionSpace>, LayeredSubdomain<FunctionSpace>>
            BC_setup(space);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(SingleSedimentaryLayer_Reg);

    static void SingleSedimentaryLayer3D(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Comm = utopia::PetscCommunicator;
        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType = FunctionSpace::SizeType;
        // using ProblemType = utopia::TensionSplitGenericPhaseField<FunctionSpace, Dim, AT1>;
        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;
        // using ProblemType = utopia::IsotropicGenericPhaseField<FunctionSpace, Dim, AT1>;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        stats.stop_and_collect("space-creation");

        if (mpi_world_rank() == 0) std::cout << "Starting Single Layer Model" << std::endl;

        stats.start();

        DamagedSedimentaryLayers<FunctionSpace> IC_setup(space, 0.0);
        DirichletAndVolConstraints<FunctionSpace, SedimentaryLayers_BC<FunctionSpace>, LayeredSubdomain<FunctionSpace>>
            BC_setup(space);

        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

        time_stepper.read(in);
        time_stepper.run();

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(SingleSedimentaryLayer3D);

//    static void DoubleSedimentaryLayer(Input &in) {
//        static const int Dim = 2;
//        static const int NVars = Dim + 1;

//        using Comm = utopia::PetscCommunicator;
//        using Mesh = utopia::PetscStructuredGrid<Dim>;
//        using Elem = utopia::PetscUniformQuad4;
//        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
//        // using SizeType = FunctionSpace::SizeType;
//        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;
//        //        using ProblemType = utopia::VolDevGenericPhaseField<FunctionSpace, Dim, AT1>;

//        Comm world;

//        MPITimeStatistics stats(world);
//        stats.start();

//        FunctionSpace space;
//        space.read(in);
//        stats.stop_and_collect("space-creation");

//        if (mpi_world_rank() == 0) std::cout << "Starting Hobbs Model" << std::endl;

//        stats.start();

//        DamagedSedimentaryLayers<FunctionSpace> IC_setup(space, 0.0);
//        DirichletAndVolConstraints<FunctionSpace, SedimentaryLayers_BC<FunctionSpace>, LayeredSubdomain<FunctionSpace>>
//            BC_setup(space);

//        IncrementalLoading<FunctionSpace, ProblemType> time_stepper(space, IC_setup, BC_setup);

//        time_stepper.read(in);
//        time_stepper.run();

//        stats.stop_collect_and_restart("end");

//        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
//        stats.stop_and_collect("output");
//        stats.describe(std::cout);
//    }

//    UTOPIA_REGISTER_APP(DoubleSedimentaryLayer);

}  // namespace utopia

#endif
