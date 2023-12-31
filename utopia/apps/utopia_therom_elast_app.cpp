
#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_WeakThermoElasticity.hpp"

#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_impl.hpp"

#include <cmath>

namespace utopia {

    template <class FunctionSpace>
    static void temperature_BC(FunctionSpace &space, const int temperature) {
        using Point = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;

        space.emplace_dirichlet_condition(
            SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 50.; }, temperature);

        for (auto s : space.mesh().sides()) {
            if (s != SideSet::bottom()) {
                space.emplace_dirichlet_condition(
                    s, UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, temperature);
            }
        }
    }

    template <class FunctionSpace>
    static void displacement_BC(FunctionSpace &space, const int disp_begin, const int disp_end) {
        using Point = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;

        for (int i = disp_begin; i < disp_end; ++i) {
            space.emplace_dirichlet_condition(
                SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.; }, i);

            space.emplace_dirichlet_condition(
                SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.; }, i);
        }
    }

    template <class FunctionSpace>
    static void thermo_elast(FunctionSpace &space, Input &in) {
        static const int Dim = FunctionSpace::Dim;
        // static const int NVars = FunctionSpace::Dim + 1;

        // expose inner types
        // using Elem = typename FunctionSpace::Shape;
        // using SizeType = typename FunctionSpace::SizeType;
        // using Dev = typename FunctionSpace::Device;
        // using Point = typename FunctionSpace::Point;

        // using Subspace = typename FunctionSpace::template Subspace<1>;
        // using ElemViewScalar = typename Subspace::ViewDevice::Elem;

        // static const int NNodes = Elem::NNodes;
        using Model = utopia::WeakThermoElasticity<FunctionSpace>;

        MPITimeStatistics stats(space.comm());

        ////////////////////////////////////////////////////////////////
        stats.start();
        // boundary conditions

        // components
        const int temperature = 0;
        const int disp_begin = 1;
        const int disp_end = disp_begin + Dim;

        bool matrix_free = true;
        in.get("matrix_free", matrix_free);

        // bool apply_BC = true;
        // in.get("apply_BC", apply_BC);

        // if(apply_BC) {
        temperature_BC(space, temperature);
        displacement_BC(space, disp_begin, disp_end);
        // }

        stats.stop_collect_and_restart("BC");

        ////////////////////////////////////////////////////////////////

        Model model(space);
        model.read(in);

        PetscMatrix H;
        PetscVector x, g;

        space.create_vector(x);
        space.create_vector(g);
        space.apply_constraints(g);

        stats.stop_collect_and_restart("model-init");

        if (matrix_free) {
            // BiCGStab<PetscMatrix, PetscVector, HOMEMADE> linear_solver;
            ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> linear_solver;
            linear_solver.verbose(true);
            linear_solver.max_it(space.n_dofs());
            linear_solver.solve(model, g, x);
            stats.stop_collect_and_restart("matrix-free-solve");

        } else {
            model.hessian(x, H);
            stats.stop_collect_and_restart("assemble");
            Factorization<PetscMatrix, PetscVector> linear_solver;
            linear_solver.solve(H, g, x);
            stats.stop_collect_and_restart("solve");
        }

        std::string output_path = "thermo_elast.vtr";

        in.get("output-path", output_path);

        rename("X", x);
        space.subspace(0).write(output_path, x);

        rename("G", g);
        space.subspace(0).write("G.vtr", g);

        stats.stop_and_collect("output");

        space.comm().root_print("n_dofs: " + std::to_string(space.n_dofs()));
        stats.describe(std::cout);
    }

    static void thermo_elast_2(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;

        FunctionSpace space;
        space.read(in);

        thermo_elast(space, in);
    }

    UTOPIA_REGISTER_APP(thermo_elast_2);

    static void thermo_elast_3(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Mesh = utopia::PetscStructuredGrid<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = FunctionSpace::SizeType;

        FunctionSpace space;
        space.read(in);

        thermo_elast(space, in);
    }

    UTOPIA_REGISTER_APP(thermo_elast_3);

}  // namespace utopia
