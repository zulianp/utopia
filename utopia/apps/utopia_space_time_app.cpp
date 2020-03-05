
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
#include "utopia_STHeatEquation.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"

#include "utopia_app_utils.hpp"

#include <cmath>

namespace utopia {

    template<class FunctionSpace, class RHSFun>
    static void space_time(
        FunctionSpace &space,
        RHSFun rhs_fun,
        Input &in)
    {
        static const int Dim   = FunctionSpace::Dim;
        static const int NVars = 1;

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
        using Model          = utopia::STHeatEquation<FunctionSpace>;

        ////Check input
        MPITimeStatistics stats(space.comm());

        ////////////////////////////////////////////////////////////////
        stats.start();
        //boundary conditions

        using Point  = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;

        ///example from: Space-time Finite Element Methods for Parabolic Initial-Boundary Problems with Variable Coefficients (Andreas Schafelner)
        bool write_mat = false;
        in.get("write_mat", write_mat);

        if(write_mat) {
            Model model(space);
            PetscVector x;
            PetscMatrix A;

            space.create_vector(x);
            space.create_matrix(A);
            model.hessian(x, A);
            rename("a", A);
            write("A.m", A);
            return;
        }

        ////////////////////////////////////////////////////////////////
        bool use_mg = false;

        if(use_mg)  {
            geometric_multigrid<Model>(space, in);
        } else {

            Model model(space);

            PetscMatrix H;
            PetscVector x, g;

            space.create_vector(x);
            space.create_vector(g);

            model.space_time_linear_form(rhs_fun, g);

            space.apply_constraints(g);
            stats.stop_collect_and_restart("model-init");

            model.hessian(x, H);
            stats.stop_collect_and_restart("assemble");

            // KSPSolver<PetscMatrix, PetscVector> linear_solver;
            // linear_solver.max_it(nx*ny);
            // linear_solver.read(in);
            Factorization<PetscMatrix, PetscVector> linear_solver;
            linear_solver.solve(H, g, x);


            stats.stop_collect_and_restart("solve");

            std::string output_path = "space_time.vtr";

            in.get("output-path", output_path);

            rename("X", x);
            space.subspace(0).write(output_path, x);


            rename("G", g);
            space.subspace(0).write("G.vtr", g);

            stats.stop_and_collect("output");

            space.comm().root_print( "n_dofs: " + std::to_string(space.n_dofs()) );
            stats.describe(std::cout);
        }
    }

    static void space_time_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = 1;

        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;
        using Scalar         = FunctionSpace::Scalar;
        using Point          = FunctionSpace::Point;

        FunctionSpace space;
        space.read(in);

        SizeType nx = 10;
        SizeType ny = 10;

        in.get("nx", nx);
        in.get("ny", ny);

        if(ny < 2.0 * (nx-1)*(nx-1)) {
            space.comm().root_print("[Warning] for avoiding oscillations discretization has to satisfy ny < 2.0 * (nx-1)*(nx-1), "
                                    "respectively dt < h^2/2");
        }

        //bottom == before / top == after
        space.emplace_dirichlet_condition(
            SideSet::bottom(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space_time(
            space,
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                static const Scalar pi = device::pi<Scalar>();
                const Scalar x = p[0];
                const Scalar y = p[1];
                return pi * device::sin(pi*x) * (device::cos(pi*y) + pi * device::sin(pi*y));
            },
            in
        );
    }

    UTOPIA_REGISTER_APP(space_time_2);

    static void space_time_3(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = 1;

        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;
        using Scalar         = FunctionSpace::Scalar;
        using Point          = FunctionSpace::Point;

        FunctionSpace space;
        space.read(in);

        SizeType nx = 10;
        SizeType ny = 10;
        SizeType nz = 10;

        in.get("nx", nx);
        in.get("ny", ny);
        in.get("nz", nz);

        if(nz < 2.0 * (nx-1)*(nx-1) || nz < 2.0 * (ny-1)*(ny-1)) {
            space.comm().root_print("[Warning] for avoiding oscillations discretization has to satisfy ny < 2.0 * (nx-1)*(nx-1), "
                                    "respectively dt < h^2/2");
        }

        //front == before / back == after
        space.emplace_dirichlet_condition(
            SideSet::front(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 10.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::bottom(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::top(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            }
        );

        space_time(
            space,
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            },
            in
        );
    }

    UTOPIA_REGISTER_APP(space_time_3);



}

