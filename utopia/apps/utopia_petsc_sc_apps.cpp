
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
//include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Core.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"

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

        using Quadrature = utopia::Quadrature<Elem, 2>;

        using SizeType = Mesh::SizeType;
        SizeType nx = 4;
        SizeType ny = 5;

        PetscCommunicator world;
        Mesh mesh(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        Quadrature quadrature;
        FunctionSpace space(mesh);

        auto &&space_view = space.view_device();
        auto &&q_view     = quadrature.view_device();

        PetscMatrix mat;
        space.create_matrix(mat);

        PhysicalGradient<FunctionSpace, Quadrature> gradients(space, quadrature);

        //View extraction
        auto &&g_view = gradients.view_device();

        Differential<FunctionSpace, Quadrature> differentials(space, quadrature);
        auto &&d_view = differentials.view_device();

        Chrono c;
        c.start();

        std::stringstream ss;

        {
            auto mat_view = device_view(mat);

            Dev::parallel_for(
                space.local_element_range(),
                [&](const SizeType &i)
            {
                Elem e;
                DofIndex dofs;
                ElementMatrix el_mat;

                space_view.elem(i, e);

                //element-wise extraction
                const auto grad = g_view.make(i, e);
                const auto dx   = d_view.make(i, e);

                el_mat.set(0.0);

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
                const SizeType n_dofs = dofs.size();

                for(SizeType i = 0; i < n_dofs; ++i) {
                    for(SizeType j = 0; j < n_dofs; ++j) {
                        mat_view.atomic_add(dofs[i], dofs[j], el_mat(i, j));
                    }
                }

                // ss << i << ")\n";
                // mesh.nodes_local(i, lnodes);
                // mesh.nodes(i, gnodes);

                // for(SizeType k = 0; k < lnodes.size(); ++k) {
                //     ss <<  "[" << lnodes[k] << ", " << gnodes[k] << " " << mesh.is_ghost(gnodes[k]) << "] ";
                // }

                // ss << "\nt = ";
                // e.describe(ss);

                // ss << std::endl;
            });
        }

        c.stop();
        std::cout << " assemblies " << c << std::endl;

        disp(mat);

        // int size = world.size();
        // int rank = world.rank();

        // world.barrier();

        // for(int i = 0; i < size; ++i) {
        //     if(i == rank) {
        //         std::cout << "--------------------------------------------\n";
        //         std::cout << ss.str() << std::endl;
        //         std::cout << "--------------------------------------------\n";
        //     }

        //     world.barrier();
        // }
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble);
}

#endif //WITH_TRILINOS
