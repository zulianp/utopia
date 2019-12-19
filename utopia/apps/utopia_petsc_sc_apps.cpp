
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
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

#include <cmath>

namespace utopia {
    static void petsc_bratu()
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;
        using Scalar           = Mesh::Scalar;
        using Point            = Mesh::Point;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 10;
        SizeType ny = scale * 10;
        SizeType nz = scale * 10;

        Mesh mesh(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        FunctionSpace space(mesh);

        //boundary conditions
        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            }
        );

        ///////////////////////////////////////

        // BratuFE<FunctionSpace> fe_model(space);
        PoissonFE<FunctionSpace> fe_model(space);
        fe_model.init_forcing_function(UTOPIA_LAMBDA(const Point &p) {
            return device::exp(-10.0 * device::abs(p[1] - 0.5));
        });

        PetscVector x;
        space.create_vector(x);

        Newton<PetscMatrix, PetscVector> newton;
        newton.verbose(true);

        space.apply_constraints(x);
        newton.solve(fe_model, x);

        rename("x", x);
        space.write("fe.vtk", x);
    }

    UTOPIA_REGISTER_APP(petsc_bratu);

    static void petsc_dm_app()
    {
        using Mesh = utopia::PetscDM<2>;
        using SizeType = Mesh::SizeType;
        SizeType nx = 10;
        SizeType ny = 10;

        PetscCommunicator world;
        Mesh dm(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        PetscMatrix mat;
        dm.create_matrix(mat);

        dm.each_element([](const Mesh::Elem &e) {
            // std::cout << e.idx() << std::endl;
        });

        std::stringstream ss;

        dm.each_node([&ss](const Mesh::Node &node) {
            assert(!node.is_ghost());
        });

        dm.each_node_with_ghosts([&ss](const Mesh::Node &node) {
            ss << "(" << node.idx() <<  ", " << node.is_ghost() << ") ";
        });

        int size = world.size();
        int rank = world.rank();

        world.barrier();

        for(int i = 0; i < size; ++i) {
            if(i == rank) {
                std::cout << "--------------------------------------------\n";
                std::cout << ss.str() << std::endl;
                std::cout << "--------------------------------------------\n";
            }

            world.barrier();
        }
    }

    UTOPIA_REGISTER_APP(petsc_dm_app);

    template<class FunctionSpace>
    static void poisson_problem(FunctionSpace &space, const bool use_direct_solver)
    {
        using Mesh             = typename FunctionSpace::Mesh;
        using Elem             = typename FunctionSpace::Elem;
        using Dev              = typename FunctionSpace::Device;
        using Vector           = typename FunctionSpace::Vector;
        using Matrix           = typename FunctionSpace::Matrix;
        using Comm             = typename FunctionSpace::Comm;

        static const int Dim    = Elem::Dim;
        static const int NNodes = Elem::NNodes;

        using DevFunctionSpace = typename FunctionSpace::ViewDevice;
        using DofIndex         = typename DevFunctionSpace::DofIndex;
        using Point            = typename FunctionSpace::Point;
        using Scalar           = typename FunctionSpace::Scalar;
        using SizeType         = typename FunctionSpace::SizeType;
        using Quadrature       = utopia::Quadrature<Elem, 2>;
        using ElementMatrix    = utopia::StaticMatrix<Scalar, NNodes, NNodes>;
        using ElementVector    = utopia::StaticVector<Scalar, NNodes>;
        using DirichletBC      = utopia::DirichletBoundaryCondition<FunctionSpace>;

        Comm &comm = space.comm();

        MPITimeStatistics stats(comm);

        Quadrature quadrature;
        auto &&space_view = space.view_device();

        comm.barrier();
        stats.start();

        Matrix mat, mass_mat;
        space.create_matrix(mat);
        space.create_matrix(mass_mat);

        Vector rhs;
        space.create_vector(rhs);

        stats.stop_and_collect("create-matrix");

        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 1.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return -1.0;
            }
        );

        auto diffusivity = UTOPIA_LAMBDA(const Point &p) -> Scalar {
            Scalar dist = 0.0;
            for(SizeType i = 0; i < p.size(); ++i) {
                Scalar v = p[i] - 0.5;
                dist += v*v;
            }

            if(device::sqrt(dist) > 0.2) {
                return 1.0;
            } else {
                return 1e-4;
            }
        };

        stats.start();


        Laplacian<FunctionSpace, Quadrature> laplacian(space, quadrature);
        MassMatrix<FunctionSpace, Quadrature> mass_matrix(space, quadrature);

        {
            //GPU assembly mock-prototype

            auto mat_view      = space.assembly_view_device(mat);
            auto mass_mat_view = space.assembly_view_device(mass_mat);
            auto rhs_view      = space.assembly_view_device(rhs);

            auto l_view = laplacian.view_device();
            auto m_view = mass_matrix.view_device();

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
                l_view.add(i, e, el_mat);
                e.centroid(c);
                el_mat *= diffusivity(c);
                space_view.add_matrix(e, el_mat, mat_view);

                //Assemble local mass-matrix and reuse el_mat
                el_mat.set(0.0);
                m_view.add(i, e, el_mat);
                space_view.add_matrix(e, el_mat, mass_mat_view);
            });
        }

        stats.stop_and_collect("assemblies");

        stats.start();
        rhs = mass_mat * rhs;

        Scalar vol = sum(mass_mat);
        std::cout << "vol: " << vol << std::endl;

        Scalar zero = sum(mat);
        std::cout << "zero: " << zero << std::endl;

        space.apply_constraints(mat, rhs);

        stats.stop_and_collect("boundary conditions ");

        stats.start();
        Vector x = rhs;
        x.set(0.0);

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

        stats.start();
        rename("x", x);
        space.write("X.vtk", x);
        stats.stop_and_collect("write");

        if(comm.rank() == 0) std::cout << "n_dofs: " << space.n_dofs() << std::endl;

        stats.describe(std::cout);
    }

    static void petsc_dm_assemble_2()
    {
        static const int Dim = 2;
        static const int NNodes = 4;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 10;
        SizeType ny = scale * 10;
        SizeType nz = 10;

        Mesh mesh(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        FunctionSpace space(mesh);
        poisson_problem(space, true);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble_2);

    static void petsc_dm_assemble_3()
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 15;
        SizeType ny = scale * 15;
        SizeType nz = scale * 15;

        Mesh mesh(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        FunctionSpace space(mesh);
        poisson_problem(space, false);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble_3);
}

#endif //WITH_TRILINOS
