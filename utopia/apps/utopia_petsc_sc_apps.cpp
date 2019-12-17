
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


    static void petsc_dm_assemble()
    {
        static const int Dim = 2;
        static const int NNodes = 4;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using Dev              = FunctionSpace::Device;

        using DevFunctionSpace = FunctionSpace::ViewDevice;
        using DofIndex         = DevFunctionSpace::DofIndex;
        using ElementMatrix    = utopia::StaticMatrix<double, NNodes, NNodes>;
        using ElementVector    = utopia::StaticVector<double, NNodes>;
        using Point            = FunctionSpace::Point;
        using Scalar           = FunctionSpace::Scalar;

        using Quadrature = utopia::Quadrature<Elem, Dim>;

        PetscCommunicator world;

        using SizeType = Mesh::SizeType;
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

        Quadrature quadrature;
        FunctionSpace space(mesh);

        auto &&space_view = space.view_device();
        auto &&q_view     = quadrature.view_device();

        world.barrier();
        c.start();

        PetscMatrix mat, mass_mat;
        space.create_matrix(mat);
        mat *= 0.0;

        space.create_matrix(mass_mat);
        mass_mat *= 0.0;

        PetscVector rhs;
        space.create_vector(rhs);
        rhs.set(1.0);

        world.barrier();
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
                    return 0.0;
                }
        ));

        bcs.push_back(
            utopia::make_unique<BoundaryCondition<FunctionSpace>>(
                space,
                SideSet::top(),
                [](const Point &p) -> Scalar {
                    return 0.0;
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
        world.barrier();

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

        rhs = mass_mat * rhs;

        for(const auto &bc : bcs) {
            bc->apply(mat, rhs);
        }

        world.barrier();
        c.stop();
        std::cout << " assemblies " << n_assemblies << " " << c << std::endl;

        SizeType nnz = utopia::nnz(mat, 0.);

        rename("r", rhs);
        space.write("R.vtk", rhs);
        std::cout << "nnz " << nnz << std::endl;


        c.start();
        PetscVector x = rhs;
        x.set(0.0);

        ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> cg;
        auto prec = std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>();
        cg.set_preconditioner(prec);
        // cg.verbose(true);
        cg.max_it(nx*ny);
        cg.rtol(1e-8);
        cg.solve(mat, rhs, x);

        c.stop();

        std::cout << c << std::endl;

        rename("x", x);
        space.write("X.vtk", x);


        // x.set(0.0);
        // for(const auto &bc : bcs) {
        //     bc->set_boundary_id(x);
        // }

        // rename("b", x);
        // space.write("B.vtk", x);


        // x.set(world.rank());

        // rename("c", x);
        // space.write("C.vtk", x);


        // each_write(x, [](const SizeType &i) {
        //     return i;
        // });

        // rename("n", x);
        // space.write("N.vtk", x);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble);
}

#endif //WITH_TRILINOS
