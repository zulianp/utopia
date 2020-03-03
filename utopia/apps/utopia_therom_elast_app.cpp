
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
#include "utopia_WeakThermoElasticity.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"

#include <cmath>

namespace utopia {

    template<class FunctionSpace>
    static void temperature_BC(
        FunctionSpace &space,
        const int temperature,
        Input &in
    )
    {
        using Point  = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;

        space.emplace_dirichlet_condition(
            SideSet::bottom(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 50.;
            },
            temperature
        );

        for(auto s : space.mesh().sides()) {
            if(s != SideSet::bottom()) {

                space.emplace_dirichlet_condition(
                    s,
                    UTOPIA_LAMBDA(const Point &) -> Scalar {
                        return 0.0;
                    },
                    temperature
                );

            }
        }
    }


    template<class FunctionSpace>
    static void displacement_BC(
        FunctionSpace &space,
        const int disp_begin,
        const int disp_end,
        Input &in
    )
    {
        using Point  = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;

        for(int i = disp_begin; i < disp_end; ++i) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return 0.;
                },
                i
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return 0.;
                },
                i
            );
        }
    }


    template<class FunctionSpace>
    static void therom_elast(
        FunctionSpace &space,
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
        using Model          = utopia::WeakThermoElasticity<FunctionSpace>;


        MPITimeStatistics stats(space.comm());

        ////////////////////////////////////////////////////////////////
        stats.start();
        //boundary conditions

        //components
        const int temperature = 0;
        const int disp_begin = 1;
        const int disp_end   = disp_begin + Dim;


        bool apply_BC = true;
        in.get("apply_BC", apply_BC);

        if(apply_BC) {
            temperature_BC(space, temperature, in);
            displacement_BC(space, disp_begin, disp_end, in);
        }

        stats.stop_collect_and_restart("BC");

        ////////////////////////////////////////////////////////////////

        Model model(space);
        model.read(in);

        PetscMatrix H;
        PetscVector x, g;

        space.create_vector(x);
        space.create_vector(g);

        stats.stop_collect_and_restart("model-init");
        model.hessian(x, H);
        // model.gradient(x, g);

        // rename("h", H);
        // write("H.m", H);


        if(!apply_BC) return;

        stats.stop_collect_and_restart("assemble");

        space.apply_constraints(g);

        Factorization<PetscMatrix, PetscVector> linear_solver;
        linear_solver.solve(H, g, x);

        stats.stop_collect_and_restart("solve");

        std::string output_path = "thermo_elast.vtr";

        in.get("output-path", output_path);

        rename("X", x);
        space.subspace(0).write(output_path, x);


        rename("G", g);
        space.subspace(0).write("G.vtr", g);

        stats.stop_and_collect("output");

        space.comm().root_print( "n_dofs: " + std::to_string(space.n_dofs()) );
        stats.describe(std::cout);
    }

    static void thermo_elast_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        FunctionSpace space;
        space.read(in);

        therom_elast(
            space,
            in
        );
    }

    UTOPIA_REGISTER_APP(thermo_elast_2);

    static void thermo_elast_3(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        FunctionSpace space;
        space.read(in);

        therom_elast(
            space,
            in
        );
    }

    UTOPIA_REGISTER_APP(thermo_elast_3);

}

