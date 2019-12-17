
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
        using Mesh             = utopia::PetscDM<2>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using Dev              = FunctionSpace::Device;

        using DevFunctionSpace = FunctionSpace::ViewDevice;
        using DofIndex         = DevFunctionSpace::DofIndex;
        using ElementMatrix    = utopia::StaticMatrix<double, 4, 4>;
        using ElementVector    = utopia::StaticVector<double, 4>;
        using Point            = FunctionSpace::Point;
        using Scalar           = FunctionSpace::Scalar;

        using Quadrature = utopia::Quadrature<Elem, 2>;

        using SizeType = Mesh::SizeType;
        SizeType scale = 20;
        SizeType nx = scale * 5;
        SizeType ny = scale * 5;

        PetscCommunicator world;

        Chrono c;
        world.barrier();
        c.start();

        Mesh mesh(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
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

        PetscMatrix mat;
        space.create_matrix(mat);

        PetscVector rhs;
        space.create_vector(rhs);
        rhs.set(0.0);

        world.barrier();
        c.stop();
        std::cout << "create-matrix: " << c << std::endl;

        PhysicalGradient<FunctionSpace, Quadrature> gradients(space, quadrature);

        //View extraction
        auto &&g_view = gradients.view_device();

        Differential<FunctionSpace, Quadrature> differentials(space, quadrature);
        auto &&d_view = differentials.view_device();

        BoundaryCondition<FunctionSpace> bc_left(
            space,
            SideSet::LEFT,
            [](const Point &p) -> Scalar {
                return -p[1];
            }
        );

        BoundaryCondition<FunctionSpace> bc_right(
            space,
            SideSet::RIGHT,
            [](const Point &p) -> Scalar {
                return p[1];
            }
        );

        c.start();
        world.barrier();

        std::stringstream ss;
        SizeType n_assemblies = 0;

        {
            auto mat_view = device_view(mat);
            //FIXME
            Write<PetscVector> w(rhs, utopia::GLOBAL_ADD);

            Dev::parallel_for(
                space.local_element_range(),
                [&](const SizeType &i)
            {
                Elem e;
                DofIndex dofs;
                ElementMatrix el_mat;
                ElementVector el_vec;

                space_view.elem(i, e);

                //element-wise extraction
                const auto grad = g_view.make(i, e);
                const auto dx   = d_view.make(i, e);

                el_mat.set(0.0);
                el_vec.set(0.0);

                const auto n = grad.n_points();
                for(std::size_t k = 0; k < n; ++k) {
                    for(std::size_t j = 0; j < grad.n_functions(); ++j) {
                        const auto g_test = grad(j, k);
                        el_mat(j, j) += dot(g_test, g_test) * dx(k);
                        for(std::size_t l = j + 1; l < grad.n_functions(); ++l) {
                            const auto v = dot(g_test, grad(l, k)) * dx(k);
                            el_mat(j, l) += v;
                            el_mat(l, j) += v;
                        }
                    }
                }

                space_view.dofs(i, dofs);

                bc_left.apply(e,  dofs, el_mat, el_vec);
                bc_right.apply(e, dofs, el_mat, el_vec);

                const SizeType n_dofs = dofs.size();

                for(SizeType i = 0; i < n_dofs; ++i) {
                    rhs.c_add(dofs[i], el_vec(i));
                    for(SizeType j = 0; j < n_dofs; ++j) {
                        mat_view.atomic_add(dofs[i], dofs[j], el_mat(i, j));
                    }
                }

                ++n_assemblies;
            });
        }

        world.barrier();
        c.stop();
        std::cout << " assemblies " << n_assemblies << " " << c << std::endl;

        SizeType nnz = utopia::nnz(mat, 0.);

        rename("a", mat);
        write("A.m", mat);
        rename("r", rhs);
        write("R.m", rhs);
        std::cout << "nnz " << nnz << std::endl;
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble);
}

#endif //WITH_TRILINOS
