
#include "utopia_Base.hpp"
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

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 4;
        SizeType ny = scale * 4;

        in.get("nx", nx);
        in.get("ny", ny);

        FunctionSpace space;

        space.build(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
            );

        space.mesh().set_field_name(0, "c");
        space.mesh().set_field_name(1, "disp_x");
        space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        stats.start();

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);  
        PFFracTension2D<FunctionSpace> BC_setup(space); 
        IncrementalLoading<FunctionSpace > time_stepper(space, IC_setup, BC_setup); 

        time_stepper.template run<IsotropicPhaseFieldForBrittleFractures<FunctionSpace>>(in); 

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

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 4;
        SizeType ny = scale * 4;

        in.get("nx", nx);
        in.get("ny", ny);

        FunctionSpace space;

        space.build(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
            );

        space.mesh().set_field_name(0, "c");
        space.mesh().set_field_name(1, "disp_x");
        space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        stats.start();

        InitialCondidtionPFTbar<FunctionSpace> IC_setup(space, 0.0);  
        PFFracFixAllDisp2D<FunctionSpace> BC_setup(space); 
        IncrementalLoading<FunctionSpace > time_stepper(space, IC_setup, BC_setup); 

        time_stepper.template run<IsotropicPhaseFieldForBrittleFractures<FunctionSpace>>(in); 

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    UTOPIA_REGISTER_APP(petsc_pressure_Tbar_isotropic_phase_field_2);







    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
    template<class FunctionSpace>
    static void isotropic_phase_field_fracture_sim(
        FunctionSpace &space,
        MPITimeStatistics &stats,
        Input &in)
    {
        static const int Dim   = FunctionSpace::Dim;
        static const int NVars = FunctionSpace::Dim + 1;

        //expose inner types
        using Comm           = typename FunctionSpace::Comm;
        using Mesh           = typename FunctionSpace::Mesh;
        using Elem           = typename FunctionSpace::Shape;
        using ElemView       = typename FunctionSpace::ViewDevice::Elem;
        using SizeType       = typename FunctionSpace::SizeType;
        using Scalar         = typename FunctionSpace::Scalar;
        using Dev            = typename FunctionSpace::Device;
        using Point          = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;

        static const int NNodes = Elem::NNodes;

        using FEFunction     = utopia::FEFunction<FunctionSpace>;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Parameters     = typename IsotropicPhaseFieldForBrittleFractures<FunctionSpace>::Parameters;

        stats.start();

        InitialCondidtionPFTension<FunctionSpace> IC_setup(space, 0.0);  
        PFFracTension2D<FunctionSpace> BC_setup(space); 
        IncrementalLoading<FunctionSpace > time_stepper(space, IC_setup, BC_setup); 

        time_stepper.template run<IsotropicPhaseFieldForBrittleFractures<FunctionSpace>>(in); 

        stats.stop_collect_and_restart("end");

        space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
    static void petsc_tension_isotropic_phase_field_3(Input &in)
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

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 4;
        SizeType ny = scale * 4;
        SizeType nz = scale * 4;

        in.get("nx", nx);
        in.get("ny", ny);
        in.get("nz", nz);

        FunctionSpace space;

        space.build(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
            );

        space.mesh().set_field_name(0, "c");
        space.mesh().set_field_name(1, "disp_x");
        space.mesh().set_field_name(2, "disp_y");
        space.mesh().set_field_name(3, "disp_z");

        stats.stop_and_collect("space-creation");

        isotropic_phase_field_fracture_sim(
            space,
            stats,
            in);
    }

    UTOPIA_REGISTER_APP(petsc_tension_isotropic_phase_field_3);
}

