#include "utopia_Base.hpp"

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
#include "utopia_L2Norm.hpp"

#include <cmath>

namespace utopia {

    void disc_test(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = 1;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using Quadrature       = utopia::Quadrature<Elem, 2>;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using ElemView         = FunctionSpace::ViewDevice::Elem;
        using Device           = FunctionSpace::Device;
        using Point            = Mesh::Point;
        using Scalar           = Mesh::Scalar;
        using Comm             = FunctionSpace::Comm;

        Comm comm;

        MPITimeStatistics stats(comm);

        stats.start();

        FunctionSpace space;
        space.read(in);



        //f(w) = 4*x^3 + y^2 + 3
        //lapl f(w) = 4*6*x + 2

        auto oracle = UTOPIA_LAMBDA(const Point &x) -> Scalar {
            return 4 * x[0] * x[0] * x[0] + x[1] * x[1] + 3;
        };

        auto lapl_oracle = UTOPIA_LAMBDA(const Point &x) -> Scalar {
            return (4.0 * 6.0) * x[0] + 2.0;
        };

        for(auto s : space.mesh().sides()) {
            space.emplace_dirichlet_condition(
                s,
                oracle,
                0
            );
        }

        stats.stop_and_collect("set-up");

        stats.start();

        PoissonFE<FunctionSpace> problem(space);
        problem.init_forcing_function(lapl_oracle);
        PetscVector x, g;
        PetscMatrix H;

        space.create_vector(x);

        problem.update(x);

        problem.gradient(x, g);
        problem.hessian(x, H);

        space.apply_constraints(H, g);

        stats.stop_and_collect("problem-assembly");
        stats.start();

        // Factorization<PetscMatrix, PetscVector> solver;
        KSPSolver<PetscMatrix, PetscVector> solver;
        solver.solve(H, g, x);

        stats.stop_and_collect("solve");

        stats.start();

        rename("x", x);
        space.write("X.vts", x);

        stats.stop_and_collect("io");


        stats.start();

        Scalar err = 0.0;

        {
            Quadrature q;
            PhysicalPoint<FunctionSpace, Quadrature> point(space, q);
            Differential<FunctionSpace, Quadrature> differential(space, q);

            FEFunction<FunctionSpace> x_fun(space, x);

            auto x_val  = x_fun.value(q);

            //view
            auto space_view = space.view_device();
            auto p_view     = point.view_device();
            auto x_view     = x_val.view_device();
            auto dx_view    = differential.view_device();

            Device::parallel_reduce(
                space.local_element_range(),
                UTOPIA_LAMBDA(const SizeType &i) -> Scalar
                {
                    StaticVector<Scalar, Quadrature::NPoints> c;
                    Point p;
                    ElemView e;

                    space_view.elem(i, e);

                    auto p_e = p_view.make(e);
                    auto c_e = x_view.make(e);
                    auto dx  = dx_view.make(e);

                    x_view.get(e, c);

                    Scalar el_err = 0.0;

                    for(SizeType qp = 0; qp < Quadrature::NPoints; ++qp) {
                        p_e.get(qp, p);
                        const Scalar diff = (c(qp) - oracle(p));
                        el_err += diff * diff * dx(qp);
                    }

                    return el_err;
                },
                err
            );
        }

        err = std::sqrt(space.comm().sum(err));

        stats.stop_and_collect("error-computation");
        std::cout << "n_dofs="   << x.size() << std::endl;
        std::cout << "l2_error=" << err      << std::endl;

        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(disc_test);

}