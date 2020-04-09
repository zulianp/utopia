

#include "utopia_Base.hpp"

//FIXME: this file causes nvcc to fail
#ifndef KOKKOS_ENABLE_CUDA

#include "utopia_RangeDevice.hpp"

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
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_TrustRegionVariableBound.hpp"
#include "utopia_QuasiNewtonBound.hpp"
#include "utopia_Backtracking.hpp"
#include "utopia_LBFGS.hpp"
#include "utopia_QuasiTrustRegionVariableBound.hpp"
#include "utopia_InitialCondition.hpp"
#include "utopia_IncrementalLoading.hpp"
#include "utopia_MLIncrementalLoading.hpp"
#include "utopia_BCSetup.hpp"

#include "utopia_ProjectedGaussSeidelNew.hpp"


#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#endif //WITH_PETSC

#include <random>
#include <cmath>
#include <chrono>

namespace utopia {


    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_tension_isotropic_phase_field_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);

        // space.mesh().set_field_name(0, "c");
        // space.mesh().set_field_name(1, "disp_x");
        // space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        stats.start();

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);
        PFFracTension2D<FunctionSpace> BC_setup(space);
        IncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace> > time_stepper(space, IC_setup, BC_setup);

        time_stepper.run(in);


        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    UTOPIA_REGISTER_APP(petsc_tension_isotropic_phase_field_2);




    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_pressure_Tbar_isotropic_phase_field_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);

        // space.mesh().set_field_name(0, "c");
        // space.mesh().set_field_name(1, "disp_x");
        // space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        stats.start();

        InitialCondidtionPFTbar<FunctionSpace> IC_setup(space, 0.0);
        PFFracFixAllDisp<FunctionSpace> BC_setup(space);
        IncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace> > time_stepper(space, IC_setup, BC_setup);

        time_stepper.run(in);


        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    UTOPIA_REGISTER_APP(petsc_pressure_Tbar_isotropic_phase_field_2);



    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_pressure_Tbar_isotropic_phase_field_2_rmtr(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);

        // space.mesh().set_field_name(0, "c");
        // space.mesh().set_field_name(1, "disp_x");
        // space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        stats.start();

        // InitialCondidtionPFTbar<FunctionSpace> IC_setup(space, 0.0);
        // PFFracFixAllDisp<FunctionSpace> BC_setup(space);

        // const auto n_levels = 3;
        // MLIncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
        //                     PFFracFixAllDisp<FunctionSpace>, InitialCondidtionPFTbar<FunctionSpace> > time_stepper(space, n_levels);


        auto n_levels = 2;
        in.get("n_levels", n_levels);

       MLIncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                            PFFracFixAllDisp<FunctionSpace>, InitialCondidtionPFFracNet<FunctionSpace> > time_stepper(space, n_levels);



        time_stepper.run(in);


        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    UTOPIA_REGISTER_APP(petsc_pressure_Tbar_isotropic_phase_field_2_rmtr);







    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_tension_phase_field_2_rmtr(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);
        // space.mesh().set_field_name(0, "c");
        // space.mesh().set_field_name(1, "disp_x");
        // space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        stats.start();

        auto n_levels = 2;
        in.get("n_levels", n_levels);


        MLIncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                            PFFracTension2D<FunctionSpace>, InitialCondidtionPFTension<FunctionSpace> > time_stepper(space, n_levels);





        time_stepper.run(in);


        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    UTOPIA_REGISTER_APP(petsc_tension_phase_field_2_rmtr);




    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_pressure_network_isotropic_phase_field_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;
        using Vector         = typename FunctionSpace::Vector;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);

        // space.mesh().set_field_name(0, "c");
        // space.mesh().set_field_name(1, "disp_x");
        // space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");
        stats.start();

        InitialCondidtionPFFracNet<FunctionSpace> IC_setup(space, 0.0);
        PFFracFixAllDisp<FunctionSpace> BC_setup(space);

        IncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace> > time_stepper(space, IC_setup, BC_setup);
        time_stepper.run(in);

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_pressure_network_isotropic_phase_field_2);


    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    static void petsc_isotropic_phase_field_3d(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        FunctionSpace space;
        space.read(in);


        // space.mesh().set_field_name(0, "c");
        // space.mesh().set_field_name(1, "disp_x");
        // space.mesh().set_field_name(2, "disp_y");
        // space.mesh().set_field_name(3, "disp_z");


        auto n_levels = 2;
        in.get("n_levels", n_levels);

        // MLIncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
        //             PFFracFixAllDisp<FunctionSpace>, InitialCondidtionPFSneddon<FunctionSpace> > time_stepper(space, n_levels);


        MLIncrementalLoading<FunctionSpace, IsotropicPhaseFieldForBrittleFractures<FunctionSpace>,
                    PFFracFixAllDisp<FunctionSpace>, InitialCondidtionPFFracNet3D<FunctionSpace> > time_stepper(space, n_levels);


        time_stepper.run(in);


        stats.stop_collect_and_restart("end");
        space.comm().root_print(std::to_string(space.n_dofs()) + " corase dofs");
        stats.stop_and_collect("space-creation");

    }

    UTOPIA_REGISTER_APP(petsc_isotropic_phase_field_3d);
}

#endif
