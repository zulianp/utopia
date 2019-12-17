
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

#include <cmath>

namespace utopia {

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
    static void poisson_problem(FunctionSpace &space)
    {
        using Mesh             = typename FunctionSpace::Mesh;
        using Elem             = typename FunctionSpace::Elem;
        using Dev              = typename FunctionSpace::Device;

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

        PetscCommunicator &comm = space.comm();

        Chrono c;

        Quadrature quadrature;

        auto &&space_view = space.view_device();
        auto &&q_view     = quadrature.view_device();

        comm.barrier();
        c.start();

        PetscMatrix mat, mass_mat;
        space.create_matrix(mat);
        mat *= 0.0;

        space.create_matrix(mass_mat);
        mass_mat *= 0.0;

        PetscVector rhs;
        space.create_vector(rhs);
        rhs.set(1.0);

        comm.barrier();
        c.stop();
        std::cout << "create-matrix: " << c << std::endl;

        PhysicalGradient<FunctionSpace, Quadrature> gradients(space, quadrature);

        //View extraction
        auto &&g_view = gradients.view_device();

        Differential<FunctionSpace, Quadrature> differentials(space, quadrature);
        auto &&d_view = differentials.view_device();

        ShapeFunction<FunctionSpace, Quadrature> functions(space, quadrature);
        auto &&f_view = functions.view_device();

        std::vector<std::unique_ptr<BoundaryCondition<FunctionSpace>>> bcs;

        bcs.push_back(
            utopia::make_unique<BoundaryCondition<FunctionSpace>>(
                space,
                SideSet::bottom(),
                [](const Point &p) -> Scalar {
                    return p[0];
                }
        ));

        bcs.push_back(
            utopia::make_unique<BoundaryCondition<FunctionSpace>>(
                space,
                SideSet::top(),
                [](const Point &p) -> Scalar {
                    return p[0];
                }
        ));

        bcs.push_back(
            utopia::make_unique<BoundaryCondition<FunctionSpace>>(
                space,
                SideSet::left(),
                [](const Point &p) -> Scalar {
                    return 0.0;
                }
        ));

        bcs.push_back(
            utopia::make_unique<BoundaryCondition<FunctionSpace>>(
                space,
                SideSet::right(),
                [](const Point &p) -> Scalar {
                    return 0.0;
                }
        ));

        c.start();
        comm.barrier();

        std::stringstream ss;
        SizeType n_assemblies = 0;

        {
            auto mat_view = device_view(mat);
            auto mass_mat_view = device_view(mass_mat);
            //FIXME
            Write<PetscVector> w(rhs, utopia::GLOBAL_ADD);

            Dev::parallel_for(
                space.local_element_range(),
                [&](const SizeType &i)
            {
                Elem e;
                DofIndex dofs;
                ElementMatrix el_mat;
                ElementMatrix el_mass_mat;
                ElementVector el_vec;

                space_view.elem(i, e);

                //element-wise extraction
                const auto grad = g_view.make(i, e);
                const auto dx   = d_view.make(i, e);
                const auto f    = f_view.make(i, e);

                el_mat.set(0.0);
                el_vec.set(0.0);
                el_mass_mat.set(0.0);

                const auto n = grad.n_points();
                for(std::size_t k = 0; k < n; ++k) {
                    for(std::size_t j = 0; j < grad.n_functions(); ++j) {
                        const auto g_test = grad(j, k);
                        el_mat(j, j) += dot(g_test, g_test) * dx(k);
                        el_mass_mat(j, j) += f(j, k) * f(j, k) * dx(k);
                        for(std::size_t l = j + 1; l < grad.n_functions(); ++l) {
                            const auto v = dot(g_test, grad(l, k)) * dx(k);
                            el_mat(j, l) += v;
                            el_mat(l, j) += v;

                            const auto v_mass = f(j, k) * f(l, k) * dx(k);
                            el_mass_mat(j, j) += v_mass;
                        }
                    }
                }

                space_view.dofs(i, dofs);

                const SizeType n_dofs = dofs.size();

                for(SizeType i = 0; i < n_dofs; ++i) {
                    rhs.c_add(dofs[i], el_vec(i));
                    for(SizeType j = 0; j < n_dofs; ++j) {
                        mat_view.atomic_add(dofs[i], dofs[j], el_mat(i, j));
                        mass_mat_view.atomic_add(dofs[i], dofs[j], el_mass_mat(i, j));
                    }
                }

                ++n_assemblies;
            });
        }

        PetscVector rhs_raw = rhs;
        rhs = mass_mat * rhs;

        for(const auto &bc : bcs) {
            bc->apply(mat, rhs);
            bc->apply(mat, rhs_raw);
        }

        comm.barrier();
        c.stop();
        std::cout << " assemblies " << n_assemblies << " " << c << std::endl;

        SizeType nnz = utopia::nnz(mat, 0.);

        rename("r", rhs);
        space.write("R.vtk", rhs);
        std::cout << "nnz " << nnz << std::endl;


        c.start();
        PetscVector x = rhs;
        x.set(0.0);

        rename("a", mat);
        write("A.m", mat);
        write("R.m", rhs);

        Factorization<PetscMatrix, PetscVector> solver;
        solver.solve(mat, rhs, x);

        // ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> cg;
        // auto prec = std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>();
        // cg.set_preconditioner(prec);
        // // cg.verbose(true);
        // cg.max_it(space.n_dofs());
        // cg.rtol(1e-8);
        // cg.solve(mat, rhs, x);

        c.stop();

        std::cout << c << std::endl;

        rename("x", x);
        space.write("X.vtk", x);

        PetscVector mass_vector = sum(mass_mat, 1);
        rhs /= mass_vector;

        rename("r", rhs_raw);
        space.write("R.vtk", rhs_raw);

        x.set(0.0);
        for(const auto &bc : bcs) {
            bc->set_boundary_id(x);
        }

        rename("b", x);
        space.write("B.vtk", x);

    }


    static void petsc_dm_assemble_2()
    {
        static const int Dim = 2;
        static const int NNodes = 4;

        //Types
        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 10;
        SizeType ny = scale * 10;
        SizeType nz = 10;

        Chrono c;
        world.barrier();
        c.start();

        Mesh mesh(
            world,
            {nx, ny},
            {1.0, 0.0},
            {2.0, 1.0}
        );

        world.barrier();
        c.stop();

        std::cout << "mesh-gen: " << c << std::endl;

        FunctionSpace space(mesh);
        poisson_problem(space);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble_2);

    static void petsc_dm_assemble_3()
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        //Types
        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 5;
        SizeType ny = scale * 5;
        SizeType nz = scale * 5;

        Chrono c;
        world.barrier();
        c.start();

        Mesh mesh(
            world,
            {nx, ny, nz},
            {1.0, 0.0, 0.0},
            {2.0, 1.0, 1.0}
        );

        world.barrier();
        c.stop();

        std::cout << "mesh-gen: " << c << std::endl;

        FunctionSpace space(mesh);
        poisson_problem(space);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble_3);
}

#endif //WITH_TRILINOS
