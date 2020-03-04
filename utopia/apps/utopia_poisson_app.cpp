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
#include "utopia_app_utils.hpp"

#include <cmath>

namespace utopia {

    static void poisson_mg(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = 1;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using Point            = FunctionSpace::Point;
        using Scalar           = FunctionSpace::Scalar;

        FunctionSpace space;
        space.read(in);


        for(int c = 0; c < space.n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return p[1];
                },
                c
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return -p[1];
                },
                c
            );
        }

        geometric_multigrid<PoissonFE<FunctionSpace>>(space, in);
    }

    UTOPIA_REGISTER_APP(poisson_mg);


    template<class FunctionSpace>
    static void poisson_problem(FunctionSpace &space, Input &in)
    {
        // using Mesh             = typename FunctionSpace::Mesh;
        using Elem             = typename FunctionSpace::Elem;
        using Dev              = typename FunctionSpace::Device;
        using Vector           = typename FunctionSpace::Vector;
        using Matrix           = typename FunctionSpace::Matrix;
        using Comm             = typename FunctionSpace::Comm;

        static const int Dim    = Elem::Dim;
        static const int NFunctions = Elem::NFunctions;

        using Point            = typename FunctionSpace::Point;
        using Scalar           = typename FunctionSpace::Scalar;
        using SizeType         = typename FunctionSpace::SizeType;
        using Quadrature       = utopia::Quadrature<Elem, 2>;
        using ElementMatrix    = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;


        bool use_direct_solver = true;
        bool debug_matrices    = false;
        bool matrix_free       = false;

        in.get("use_direct_solver", use_direct_solver);
        in.get("debug_matrices",    debug_matrices);
        in.get("matrix_free",       matrix_free);


        Comm &comm = space.comm();
        MPITimeStatistics stats(comm);

        //////////////////////////////////////////
        stats.start();

        Quadrature quadrature;
        auto &&space_view = space.view_device();

        Vector rhs, x;
        space.create_vector(rhs);
        space.create_vector(x);

        stats.stop_and_collect("create-vector");
        //////////////////////////////////////////

        for(int c = 0; c < space.n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return p[1];
                },
                c
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return -p[1];
                },
                c
            );
        }

        auto diffusivity = UTOPIA_LAMBDA(const Point &p) -> Scalar {
            // Scalar dist = 0.0;
            // for(SizeType i = 0; i < p.size(); ++i) {
            //     Scalar v = p[i] - 0.5;
            //     dist += v*v;
            // }

            // if(device::sqrt(dist) > 0.2) {
            //     return 1.0;
            // } else {
            //     return 1e-4;
            // }

            return 1.0;
        };

        stats.stop_and_collect("set-up");
        //////////////////////////////////////////

        if(matrix_free) {
            //////////////////////////////////////////
            const SizeType n_iter = space.n_dofs(); assert(n_iter > 0);
            PoissonFE<FunctionSpace> poisson(space);

            space.apply_constraints(rhs);

            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            cg.verbose(true);
            cg.max_it(n_iter);
            cg.rtol(1e-8);
            cg.solve(poisson, rhs, x);

            stats.stop_and_collect("matrix-free-solve");
            //////////////////////////////////////////

        } else {

            auto lapl   = laplacian(space, quadrature);
            auto mass_m = mass_matrix(space, quadrature);

            Matrix mat, mass_mat;
            space.create_matrix(mat);
            // space.create_matrix(mass_mat);
            //copying is faster than create_matrix
            mass_mat = mat;

            stats.stop_and_collect("create-matrices");

            {
                //GPU assembly mock-prototype
                auto mat_view      = space.assembly_view_device(mat);
                auto mass_mat_view = space.assembly_view_device(mass_mat);
                auto rhs_view      = space.assembly_view_device(rhs);

                auto l_view = lapl.view_device();
                auto m_view = mass_m.view_device();

                if(debug_matrices) {
                    disp("lapl");
                    l_view.describe();

                    disp("mass_m");
                    m_view.describe();
                }

                Dev::parallel_for(
                    space.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                {
                    Elem e;

                    //FIXME this is too big for GPU stack memory for hexas
                    ElementMatrix el_mat;
                    Point c;
                    space_view.elem(i, e);

                    //Assemble local laplacian
                    el_mat.set(0.0);
                    l_view.assemble(i, e, el_mat);
                    e.centroid(c);
                    el_mat *= diffusivity(c);
                    space_view.add_matrix(e, el_mat, mat_view);

                    //Assemble local mass-matrix and reuse el_mat
                    el_mat.set(0.0);
                    m_view.assemble(i, e, el_mat);
                    space_view.add_matrix(e, el_mat, mass_mat_view);
                });
            }

            stats.stop_and_collect("assemblies");

            stats.start();
            rhs = mass_mat * rhs;

            // rename("m",  mass_mat);
            // write("M.m", mass_mat);

            // rename("a",  mat);
            // write("A.m", mat);

            Scalar vol = sum(mass_mat);
            std::cout << "vol: " << vol << std::endl;

            Scalar zero = sum(mat);
            std::cout << "zero: " << zero << std::endl;

            space.apply_constraints(mat, rhs);

            stats.stop_collect_and_restart("boundary conditions ");

            if(use_direct_solver) {
                Factorization<Matrix, Vector> solver;
                solver.solve(mat, rhs, x);
            }  else {
                ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
                auto prec = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
                cg.set_preconditioner(prec);
                cg.verbose(true);

                const SizeType n_iter = space.n_dofs();

                assert(n_iter > 0);

                cg.max_it(n_iter);
                cg.rtol(1e-8);
                cg.solve(mat, rhs, x);
            }

            stats.stop_and_collect("solve");
        }

        //////////////////////////////////////////
        stats.start();
        rename("x", x);
        space.write("X.vtk", x);
        stats.stop_and_collect("write");
        //////////////////////////////////////////

        if(comm.rank() == 0) std::cout << "n_dofs: " << space.n_dofs() << std::endl;
        stats.describe(std::cout);
    }

    static void petsc_poisson_2(Input &in)
    {
        static const int Dim = 2;
        static const int NNodes = 4;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;
        FunctionSpace space;
        space.read(in);

        poisson_problem(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_poisson_2);

    static void petsc_poisson_3(Input &in)
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;

        PetscCommunicator world;
        FunctionSpace space;
        space.read(in);

        poisson_problem(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_poisson_3);

    static void petsc_dm_mvar_poisson_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;

        FunctionSpace space;
        space.read(in);

        poisson_problem(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_dm_mvar_poisson_2);

    static void petsc_dm_mvar_poisson_3(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;

        FunctionSpace space;
        space.read(in);

        poisson_problem(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_dm_mvar_poisson_3);

}