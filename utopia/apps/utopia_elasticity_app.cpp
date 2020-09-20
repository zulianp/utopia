#include "utopia_Base.hpp"

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
#include "utopia_LinearElasticityFE.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PhaseField.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_app_utils.hpp"

// petsc
#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_FE.hpp"
#include "utopia_petsc_Matrix.hpp"

// std
#include <cmath>

namespace utopia {

    template <int Dim>
    using MeshType = utopia::PetscDMDA<StaticVector<PetscScalar, Dim>, ArrayView<PetscInt, Dim>>;

    // template<int Dim>
    // using MeshType = utopia::PetscStructuredGrid<Dim>;

    static void elast_mg_2(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh = utopia::MeshType<Dim>;
        // using Elem             = utopia::PetscUniformQuad4;
        using Elem = utopia::Tri3<PetscScalar, 2>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using Point = FunctionSpace::Point;
        using Scalar = FunctionSpace::Scalar;

        FunctionSpace space;
        space.read(in);

        for (int c = 0; c < FunctionSpace::n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.1 * Scalar(c == 1); }, c);

            space.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return -0.1 * Scalar(c == 1); }, c);
        }

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after init");

        geometric_multigrid<LinearElasticityFE<FunctionSpace>>(space, in);
    }

    UTOPIA_REGISTER_APP(elast_mg_2);

    static void elast_mg_3(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh = utopia::MeshType<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using Point = FunctionSpace::Point;
        using Scalar = FunctionSpace::Scalar;

        FunctionSpace space;
        space.read(in);

        for (int c = 0; c < FunctionSpace::n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.1 * Scalar(c == 1); }, c);

            space.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return -0.1 * Scalar(c == 1); }, c);
        }

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after init");

        geometric_multigrid<LinearElasticityFE<FunctionSpace>>(space, in);
    }

    UTOPIA_REGISTER_APP(elast_mg_3);

    template <class FunctionSpace>
    static void linear_elasticity(FunctionSpace &space, Input &in) {
        using Elem = typename FunctionSpace::Elem;
        // using Dev = typename FunctionSpace::Device;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Comm = typename FunctionSpace::Comm;

        // static const int Dim = Elem::Dim;
        // static const int NFunctions = Elem::NFunctions;

        using Point = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Quadrature = utopia::Quadrature<Elem, 2>;
        // using ElementMatrix = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;

        bool use_direct_solver = false;
        bool debug_matrices = false;
        bool matrix_free = false;
        std::string output_path = "elasticity.vtr";

        in.get("use_direct_solver", use_direct_solver);
        in.get("debug_matrices", debug_matrices);
        in.get("matrix_free", matrix_free);
        in.get("output_path", output_path);

        Comm &comm = space.comm();

        MPITimeStatistics stats(comm);
        stats.start();

        Quadrature quadrature;
        // auto &&space_view = space.view_device();

        //////////////////////////////////////////////////////////////////
        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("before-create-vector");

        Vector rhs, x;
        space.create_vector(rhs);
        space.create_vector(x);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-create-vector");
        stats.stop_collect_and_restart("create-vector");

        //////////////////////////////////////////////////////////////////

        for (int c = 0; c < space.n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.1 * Scalar(c == 1); }, c);

            space.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return -0.1 * Scalar(c == 1); }, c);
        }

        space.apply_constraints(rhs);

        stats.stop_collect_and_restart("setup");

        LinearElasticityFE<FunctionSpace> lin_elast(space);
        lin_elast.read(in);

        const SizeType n_iter = space.n_dofs();
        assert(n_iter > 0);

        if (matrix_free) {
            ConjugateGradient<Matrix, Vector, HOMEMADE> solver;
            solver.apply_gradient_descent_step(true);

            UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after solver allocation");
            solver.verbose(true);

            solver.max_it(n_iter);
            solver.rtol(1e-6);
            solver.atol(1e-6);
            solver.solve(lin_elast, rhs, x);

            stats.stop_collect_and_restart("matrix-free-solve");

        } else {
            //////////////////////////////////////////////////////////////////

            UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("before-create-matrix");

            Matrix mat;
            space.create_matrix(mat);
            stats.stop_collect_and_restart("create_matrix");

            UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-create-matrix");
            //////////////////////////////////////////////////////////////////

            lin_elast.hessian(x, mat);

            UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after assemblies");
            stats.stop_collect_and_restart("assemblies");

            //////////////////////////////////////////////////////////////////

            UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after boundary conditions");
            stats.stop_collect_and_restart("boundary conditions ");
            comm.root_print("Solving...");

            if (use_direct_solver) {
                Factorization<Matrix, Vector> solver;
                solver.solve(mat, rhs, x);
            } else if (x.size() > 1e6) {
                KSPSolver<Matrix, Vector> solver;
                solver.verbose(true);

                solver.max_it(n_iter);
                solver.rtol(1e-6);
                solver.atol(1e-6);

                solver.solve(mat, rhs, x);

            } else {
                ConjugateGradient<Matrix, Vector, HOMEMADE> solver;
                solver.apply_gradient_descent_step(true);

                auto prec = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
                solver.set_preconditioner(prec);
                solver.verbose(true);

                solver.max_it(n_iter);
                solver.rtol(1e-6);
                solver.atol(1e-6);
                solver.solve(mat, rhs, x);
            }

            UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after solve");
            stats.stop_and_collect("solve");
        }

        stats.start();
        rename("x", x);
        space.write(output_path, x);
        stats.stop_and_collect("write");

        comm.root_print("n_dofs: " + std::to_string(space.n_dofs()));
        stats.describe(std::cout);
    }

    static void petsc_elasticity_2(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh = utopia::MeshType<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = Mesh::SizeType;

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("start");

        FunctionSpace space;
        space.read(in);

        // space.describe();

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("space.read(in)");

        linear_elasticity(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_elasticity_2);

    static void petsc_elasticity_3(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh = utopia::MeshType<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = Mesh::SizeType;

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("start");

        FunctionSpace space;
        space.read(in);

        // space.mesh().set_field_names({
        //     "disp_x",
        //     "disp_y",
        //     "disp_z",
        // });

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("space.read(in)");

        linear_elasticity(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_elasticity_3);

    static void petsc_matrix_free_test(Input &in) {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh = utopia::MeshType<Dim>;
        using Elem = utopia::PetscUniformQuad4;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using SizeType = Mesh::SizeType;
        using Scalar = Mesh::Scalar;
        using Point = Mesh::Point;

        FunctionSpace space;
        space.read(in);

        for (int c = 0; c < FunctionSpace::n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &p)->Scalar { return c + 1.0 + p[0]; }, c);

            space.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &p)->Scalar { return c + 2.0 + p[1]; }, c);
        }

        LinearElasticityFE<FunctionSpace> lin_elast(space);
        lin_elast.read(in);

        PetscVector x, y;
        space.create_vector(x);
        space.create_vector(y);

        auto space_x = space.subspace(0);
        auto space_y = space.subspace(1);

        space_x.sample(
            x, UTOPIA_LAMBDA(const Point &x) { return -x[0] * x[0] * 100; });

        space_y.sample(
            x, UTOPIA_LAMBDA(const Point &x) { return x[1] * x[1] * 100; });

        space.apply_constraints(x);

        lin_elast.apply(x, y);
        PetscMatrix H;
        space.create_matrix(H);
        lin_elast.hessian(x, H);

        PetscVector y_mat = H * x;
        PetscVector diff = y - y_mat;

        // disp(diff);

        Scalar norm_diff = norm_infty(diff);

        disp(norm_diff);

        assert(norm_diff < 1e-10);

        space_x.sample(
            x, UTOPIA_LAMBDA(const Point &x) { return -x[0] * x[0] * 200; });

        space_y.sample(
            x, UTOPIA_LAMBDA(const Point &x) { return x[1] * x[1] * 200; });

        lin_elast.apply(x, y);

        y_mat = H * x;
        diff = y - y_mat;

        // disp(diff);

        norm_diff = norm_infty(diff);

        disp(norm_diff);

        assert(norm_diff < 1e-10);
    }

    UTOPIA_REGISTER_APP(petsc_matrix_free_test);

    static void petsc_strain(Input &in) {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh = utopia::MeshType<Dim>;
        using Elem = utopia::PetscUniformHex8;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;
        // using ElemView = FunctionSpace::ViewDevice::Elem;
        // using SizeType = Mesh::SizeType;
        // using Scalar = Mesh::Scalar;
        // using Quadrature = utopia::Quadrature<Elem, 2>;
        // using Dev = FunctionSpace::Device;
        // using VectorD = utopia::StaticVector<Scalar, Dim>;

        FunctionSpace space;
        space.read(in);

        PetscVector u;
        space.create_vector(u);
        u.set(0.1);

        // compute_strain_energy_splitting(space, u);
    }

    UTOPIA_REGISTER_APP(petsc_strain);

    static void test_elast_expr() {
        using Scalar = double;
        StaticMatrix<Scalar, 3, 3> strain, stress;
        strain.set(0.0);

        strain(0, 0) = 1.0;
        strain(1, 0) = strain(0, 1) = 0.5;
        strain(2, 0) = strain(0, 2) = 0.5;

        stress = 2 * 80.0 * strain + 120.0 * trace(strain) * (device::identity<Scalar>());
        disp(stress);
    }

    UTOPIA_REGISTER_APP(test_elast_expr);

}  // namespace utopia
