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

#include <cmath>


namespace utopia {

    template<class FunctionSpace>
    static void linear_elasticity(FunctionSpace &space, Input &in)
    {
        // using Mesh             = typename FunctionSpace::Mesh;
        using Elem             = typename FunctionSpace::Elem;
        using Dev              = typename FunctionSpace::Device;
        using Vector           = typename FunctionSpace::Vector;
        using Matrix           = typename FunctionSpace::Matrix;
        using Comm             = typename FunctionSpace::Comm;

        static const int Dim    = Elem::Dim;
        // static const int NNodes = Elem::NNodes;
        static const int NFunctions = Elem::NFunctions;

        // using DevFunctionSpace = typename FunctionSpace::ViewDevice;
        using Point            = typename FunctionSpace::Point;
        using Scalar           = typename FunctionSpace::Scalar;
        using SizeType         = typename FunctionSpace::SizeType;
        using Quadrature       = utopia::Quadrature<Elem, 2>;
        using ElementMatrix    = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;
        // using ElementVector    = utopia::StaticVector<Scalar, NFunctions>;


        bool use_direct_solver = false;
        bool debug_matrices = false;
        std::string output_path = "elasticity.vtr";

        in.get("use_direct_solver", use_direct_solver);
        in.get("debug_matrices", debug_matrices);
        in.get("output_path", output_path);

        Comm &comm = space.comm();

        MPITimeStatistics stats(comm);

        Quadrature quadrature;
        auto &&space_view = space.view_device();

        comm.barrier();
        stats.start();

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        Matrix mat, mass_mat;
        space.create_matrix(mat);
        // space.create_matrix(mass_mat);
        //copying is cheaper than create_matrix
        mass_mat = mat;

        Vector rhs;
        space.create_vector(rhs);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        stats.stop_and_collect("create-matrix");

        // compression_only_bc(space);

        for(int c = 0; c < space.n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return 0.1 * Scalar(c==1);
                },
                c
            );

            space.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return -0.1*Scalar(c==1);
                },
                c
            );
        }

        stats.start();

        LinearElasticity<FunctionSpace, Quadrature> elast(space, quadrature, 80.0, 120.0);
        MassMatrix<FunctionSpace, Quadrature> mass_matrix(space, quadrature);

        {
            auto mat_view      = space.assembly_view_device(mat);
            auto mass_mat_view = space.assembly_view_device(mass_mat);
            auto rhs_view      = space.assembly_view_device(rhs);

            auto elast_view = elast.view_device();
            auto m_view     = mass_matrix.view_device();

            if(debug_matrices) {
                disp("elast");
                elast_view.describe();

                disp("mass_matrix");
                m_view.describe();
            }

            Dev::parallel_for(
                space.local_element_range(),
                UTOPIA_LAMBDA(const SizeType &i)
            {
                Elem e;

                //FIXME this is too big for GPU stack memory for hexas
                ElementMatrix el_mat;
                space_view.elem(i, e);

                //Assemble local elast
                el_mat.set(0.0);
                elast_view.assemble(i, e, el_mat);
                space_view.add_matrix(e, el_mat, mat_view);

                //Assemble local mass-matrix and reuse el_mat
                el_mat.set(0.0);
                m_view.assemble(i, e, el_mat);
                space_view.add_matrix(e, el_mat, mass_mat_view);
            });
        }

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        stats.stop_collect_and_restart("assemblies");

        rhs = mass_mat * rhs;

        Scalar vol = sum(mass_mat);
        comm.root_print("vol: " + std::to_string(vol) );

        Scalar zero = sum(mat);
        comm.root_print("zero: " + std::to_string(zero) );

        space.apply_constraints(mat, rhs);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        stats.stop_collect_and_restart("boundary conditions ");

        Vector x = rhs;
        x.set(0.0);

        comm.root_print("Solving...");

        const SizeType n_iter = space.n_dofs(); assert(n_iter > 0);

        if(use_direct_solver) {
            Factorization<Matrix, Vector> solver;
            solver.solve(mat, rhs, x);
        } else if(x.size() > 1e6) {
            KSPSolver<Matrix, Vector> solver;
            solver.verbose(true);
            
            solver.max_it(n_iter);
            solver.rtol(1e-6);
            solver.atol(1e-6);

            solver.solve(mat, rhs, x);

        } else {
            ConjugateGradient<Matrix, Vector, HOMEMADE> solver;

            auto prec = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            solver.set_preconditioner(prec);
            solver.verbose(true);

            solver.max_it(n_iter);
            solver.rtol(1e-6);
            solver.atol(1e-6);
            solver.solve(mat, rhs, x);
        }

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        stats.stop_and_collect("solve");

        // stats.start();
        // compute_strain_energy_splitting(space, x);
        // stats.stop_and_collect("splitting");


        stats.start();
        rename("x", x);
        space.write(output_path, x);
        stats.stop_and_collect("write");

        // rename("rhs", rhs);
        // space.write("R.vtk", rhs);

        comm.root_print( "n_dofs: " + std::to_string(space.n_dofs()) );
        stats.describe(std::cout);
    }

    static void petsc_elasticity_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType         = Mesh::SizeType;

        FunctionSpace space;
        space.read(in);
        linear_elasticity(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_elasticity_2);


    static void petsc_elasticity_3(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType         = Mesh::SizeType;

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        FunctionSpace space;
        space.read(in);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE();

        linear_elasticity(space, in);
    }

    UTOPIA_REGISTER_APP(petsc_elasticity_3);

    static void petsc_strain()
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using ElemView       = FunctionSpace::ViewDevice::Elem;
        using SizeType       = Mesh::SizeType;
        using Scalar         = Mesh::Scalar;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Dev            = FunctionSpace::Device;
        using VectorD        = utopia::StaticVector<Scalar, Dim>;
        // using MatrixDxD        = utopia::StaticMatrix<Scalar, Dim, Dim>;



        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 2;
        SizeType ny = scale * 2;
        SizeType nz = scale * 2;

        FunctionSpace space;

        space.build(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        PetscVector u;
        space.create_vector(u);
        u.set(0.1);

        // compute_strain_energy_splitting(space, u);
    }

    UTOPIA_REGISTER_APP(petsc_strain);


    static void test_elast_expr()
    {
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

}
