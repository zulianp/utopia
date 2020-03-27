
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
#include "utopia_PhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_TrustRegionVariableBound.hpp"

#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"

#include <random>
#include <cmath>

namespace utopia {

    template<class FunctionSpace>
    static void init_phase_field(FunctionSpace &space, PetscVector &x)
    {
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh           = typename FunctionSpace::Mesh;
        using Elem           = typename FunctionSpace::Shape;
        using ElemView       = typename FunctionSpace::ViewDevice::Elem;
        using SizeType       = typename FunctionSpace::SizeType;
        using Scalar         = typename FunctionSpace::Scalar;
        using Dev            = typename FunctionSpace::Device;
        using Point          = typename FunctionSpace::Point;

        using Subspace       = typename FunctionSpace::template Subspace<1>;
        using ElemViewScalar = typename Subspace::ViewDevice::Elem;

        static const int NNodes = Elem::NNodes;

        // un-hard-code
        auto C = space.subspace(0);

        auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) -> Scalar {
            // const Scalar dist_x = 0.5 - x[0];
            Scalar f = 0.0;
            // for(int i = 1; i < Dim; ++i) {
                // auto dist_i = x[1];
                //f += device::exp(-500.0 * x[i] * x[i]);
                if(  x[0] > (0.5-space.mesh().min_spacing()) && x[0] < (0.5 + space.mesh().min_spacing())  && x[1]  < 0.5 ){
                    f = 1.0;
                    // f = 0.0;
                }
                else{
                    f = 0.0;
                }
            // }

            return f;
        });

        {
            auto C_view       = C.view_device();
            auto sampler_view = sampler.view_device();
            auto x_view       = space.assembly_view_device(x);

            Dev::parallel_for(space.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                ElemViewScalar e;
                C_view.elem(i, e);

                StaticVector<Scalar, NNodes> s;
                sampler_view.assemble(e, s);
                C_view.set_vector(e, s, x_view);
            });
        }
    }


    template<class FunctionSpace>
    static void enforce_BC_time_dependent(FunctionSpace &space, const typename FunctionSpace::Scalar & disp,  const typename FunctionSpace::Scalar & t){

        static const int Dim   = FunctionSpace::Dim;

        using Point          = typename FunctionSpace::Point;
        using Scalar         = typename FunctionSpace::Scalar;

        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            },
            1
            );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return t*disp;
            },
            1
            );

        for(int d = 2; d < Dim + 1; ++d) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
                },
                d
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return 0.0;
                },
                d
                );
            }
    }


    template<class FunctionSpace>
    static void build_irreversility_constraint(const PetscVector &x_old, PetscVector &x_new, const typename FunctionSpace::SizeType & comp)
    {
        static const int Dim   = FunctionSpace::Dim;
        using Scalar         = typename FunctionSpace::Scalar;
        using SizeType       = typename FunctionSpace::SizeType;

        {
            auto d_x_old = const_device_view(x_old);

            parallel_transform(x_new, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar
            {
                if(i%(Dim+1)==comp)
                    return d_x_old.get(i);
                else
                    return -9e15;

            });
        }

    }



    template<class FunctionSpace>
    static void phase_field_fracture_sim(
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

        using Subspace       = typename FunctionSpace::template Subspace<1>;
        using ElemViewScalar = typename Subspace::ViewDevice::Elem;

        static const int NNodes = Elem::NNodes;

        using FEFunction     = utopia::FEFunction<FunctionSpace>;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Parameters     = typename PhaseFieldForBrittleFractures<FunctionSpace>::Parameters;

        auto &mesh = space.mesh();

        Scalar disp = 0.001;

        in.get("disp", disp);

        bool with_damage = true;
        in.get("with_damage", with_damage);

        bool with_BC = true;
        in.get("with_BC", with_BC);

        stats.start();



        stats.stop_collect_and_restart("BC");

        // PhaseFieldForBrittleFractures<FunctionSpace> pp(space);
        // pp.read(in);


        PetscVector x;
        space.create_vector(x);
        x.set(0.0);

        if(with_damage) {
            init_phase_field(space, x);
        }

        stats.stop_collect_and_restart("phase-field-init");


        Scalar dt = 1e-4;
        Scalar time_=dt;
        Scalar num_ts = 100;
        std::string output_path = "phase_field";
        // print IG
        rename("X", x);

        PetscVector irreversibility_constraint = x;

        // as space gets copied, we need to instantiate PF problem every time BC changes ...
        PhaseFieldForBrittleFractures<FunctionSpace> pp(space);
        pp.read(in);


        space.write(output_path+"_"+std::to_string(0.0)+".vtr", x);
        for (auto t=0; t < num_ts; t++)
        {
            std::cout<<"Time-step: "<< t << "  \n";

            if(with_BC) {
                space.reset_bc();
                enforce_BC_time_dependent(space, disp, time_);
            }

            space.apply_constraints(x);

                                                                                        // PF component
            build_irreversility_constraint<FunctionSpace>(x, irreversibility_constraint, 0);

            //auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector>>();
            // Newton<PetscMatrix, PetscVector> solver(linear_solver);
            // in.get("solver", solver);

            // auto qp_solver = std::make_shared<utopia::MPGRP<PetscMatrix, PetscVector> >();
            auto qp_solver =  std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector> >(linear_solver);
            TrustRegionVariableBound<PetscMatrix, PetscVector> solver(qp_solver);
            auto box = make_lower_bound_constraints(make_ref(irreversibility_constraint));
            solver.set_box_constraints(box);
            in.get("solver", solver);
            solver.solve(pp, x);

            rename("X", x);
            space.write(output_path+"_"+std::to_string(time_)+".vtr", x);

            // increment time step
            time_+=dt;
        }


        stats.stop_collect_and_restart("solve+assemble");

        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    static void petsc_phase_field_2(Input &in)
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

        phase_field_fracture_sim(
            space,
            stats,
            in
        );
    }

    UTOPIA_REGISTER_APP(petsc_phase_field_2);

    static void petsc_phase_field_3(Input &in)
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

        stats.stop_and_collect("space-creation");

        phase_field_fracture_sim(
            space,
            stats,
            in);
    }

    UTOPIA_REGISTER_APP(petsc_phase_field_3);
}

