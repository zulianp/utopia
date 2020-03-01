
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

#include <cmath>

namespace utopia {

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
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;

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

        if(with_BC) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return -disp;
                },
                1
                );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return disp;
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

        stats.stop_collect_and_restart("BC");

        PhaseFieldForBrittleFractures<FunctionSpace> pp(space);
        pp.read(in);

        PetscMatrix H;
        PetscVector x, g;
        Scalar f;

        space.create_vector(x);

        x.set(0.0);

        auto C = space.subspace(0);

        if(with_damage) {

            auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) -> Scalar {
                const Scalar dist_x = 0.5 - x[0];
                Scalar f = device::exp(-500.0 * dist_x * dist_x);
                for(int i = 1; i < Dim; ++i) {
                    auto dist_i = x[1];
                    f += device::exp(-500.0 * x[i] * x[i]);
                }

                return f;
            });

            {
                auto C_view       = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view       = space.assembly_view_device(x);

                Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                    ElemViewScalar e;
                    C_view.elem(i, e);

                    StaticVector<Scalar, NNodes> s;
                    sampler_view.assemble(e, s);
                    C_view.set_vector(e, s, x_view);
                });
            }
        }

        space.apply_constraints(x);

        stats.stop_collect_and_restart("phase-field-init");

    // TrustRegion<PetscMatrix, PetscVector> solver;

        auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
        Newton<PetscMatrix, PetscVector> solver(linear_solver);
        in.get("solver", solver);

        //FIXME
        // solver.solve(pp, x);

        //REMOVE ME
        pp.hessian(x, H);
        pp.gradient(x, g);
        space.apply_constraints(g);
        linear_solver->solve(H, g, x);
        //


        stats.stop_collect_and_restart("solve+assemble");


        std::string output_path = "phase_field.vtr";

        in.get("output-path", output_path);

        rename("X", x);
        C.write(output_path, x);

        stats.stop_and_collect("output");

        stats.describe(std::cout);

    // rename("X", r);
    // C.write(output_path, r);
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

        phase_field_fracture_sim(
            space,
            stats,
            in);
    }

    UTOPIA_REGISTER_APP(petsc_phase_field_3);
}

