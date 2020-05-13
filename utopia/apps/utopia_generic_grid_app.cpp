
#include "utopia_Base.hpp"

#include "utopia_petsc_Base.hpp"

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)

#include "petscfe.h"
#include "utopia_Algorithms.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_Rename.hpp"
#include "utopia_StructuredGrid.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DMPlex.hpp"
#include "utopia_petsc_DMPlex_FunctionSpace.hpp"
// #include "utopia_petsc_FE.hpp"
#include "utopia_ui.hpp"

namespace utopia {
    using V = utopia::StaticVector<PetscScalar, 2>;
    using I = utopia::ArrayView<PetscInt, 2>;
    using VA = utopia::ArrayView<PetscScalar, 2>;
    template class StructuredGrid<V, I>;
    template class PetscDMDA<V, I>;
    template class PetscDMPlex<V, I>;

    // FIXME
    void generic_grid_test(Input &in) {
        // COMPILE-TIME BEGIN (everything is evaualted at compile time)
        constexpr StructuredGrid<V, I> g({10, 10}, {0, 0}, {10, 10}, {0, 0}, {10, 10}, VA({0.0, 0.0}), VA({1.0, 1.0}));

        static_assert(g.dim() == 2, "dim must be constexpr");
        static_assert(g.n_nodes() == 100, "n_nodes has to be computable at compile time");
        static_assert(g.n_elements() == 81, "wrong number of elements");

        constexpr bool is_b = g.is_node_on_boundary(0);
        static_assert(is_b, "node 0 must be on boundary");
        // constexpr bool is_b = g.is_node_on_boundary_local_no_ghost(0, SideSet::left());

        // for some reaons clang is able to do this at compile time but not gcc
        /*constexpr*/ PetscScalar meas = g.measure();
        // COMPILE-TIME END
        assert(device::approxeq(meas, 1.0, 1e-8));
    }

    UTOPIA_REGISTER_APP(generic_grid_test);

    void vec_grid_test(Input &in) {
        // run-time sizes real data
        std::vector<PetscInt> zeros_v = {0, 0};
        std::vector<PetscInt> dims_v = {10, 10};
        std::vector<PetscScalar> box_min_v = {0.0f, 0.0f};
        std::vector<PetscScalar> box_max_v = {1.0f, 1.0f};

        // views on data
        ArrayView<PetscInt> zeros(&zeros_v[0], zeros_v.size());
        ArrayView<PetscInt> dims(&dims_v[0], dims_v.size());
        ArrayView<PetscScalar> box_min(&box_min_v[0], box_min_v.size());
        ArrayView<PetscScalar> box_max(&box_max_v[0], box_max_v.size());

        StructuredGrid<VectorView<ArrayView<PetscScalar>>, ArrayView<PetscInt>> g(
            dims, zeros, dims, zeros, dims, box_min, box_max);

        assert(g.dim() == 2);
        assert(g.n_nodes() == 100);
        assert(g.n_elements() == 81);
        PetscScalar meas = g.measure();
        bool is_b = g.is_node_on_boundary(0);
        assert(is_b);

        assert(device::approxeq(meas, 1.0, 1e-8));
    }

    UTOPIA_REGISTER_APP(vec_grid_test);

    void dmda_test(Input &in) {
        PetscCommunicator comm;
        PetscDMDA<V, I> dmda(comm);

        dmda.read(in);

        PetscVector v;
        dmda.create_vector(v);
        v.set(1.0);

        dmda.write("prova.vtr", v);
    }

    UTOPIA_REGISTER_APP(dmda_test);

    void dmplex_test(Input &in) {
        using I = utopia::ArrayView<PetscInt, 4>;
        using Mesh = utopia::PetscDMPlex<V, I>;
        static const int NVar = 2;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NVar, utopia::Tri3<PetscReal>>;
        using Elem = FunctionSpace::Elem;
        using Quadrature = utopia::Quadrature<Elem, 2>;
        using Device = FunctionSpace::Device;
        using Point = FunctionSpace::Point;
        using Scalar = FunctionSpace::Scalar;
        using ElementMatrix = utopia::StaticMatrix<Scalar, 3 * NVar, 3 * NVar>;
        // using Laplacian = utopia::Laplacian<FunctionSpace, Quadrature>;

        PetscCommunicator comm;
        MPITimeStatistics stats(comm);

        ////////////////////////////////////////////////////////////

        stats.start();

        FunctionSpace space(comm);
        space.read(in);
        auto &dmplex = space.mesh();
        dmplex.describe();

        // for (int c = 0; c < space.n_components(); ++c) {
        //     space.emplace_dirichlet_condition(
        //         SideSet::left(), UTOPIA_LAMBDA(const Point &p)->Scalar { return p[1]; }, c);

        //     space.emplace_dirichlet_condition(
        //         SideSet::right(), UTOPIA_LAMBDA(const Point &p)->Scalar { return -p[1]; }, c);
        // }

        stats.stop_collect_and_restart("dmplex-setup");

        ////////////////////////////////////////////////////////////

        PetscVector v;
        dmplex.create_vector(v);
        v.set(1.0);

        v.comm().root_print("dofs = " + std::to_string(v.size()));

        stats.stop_collect_and_restart("create-vector");
        ////////////////////////////////////////////////////////////

        utopia::rename("X", v);
        dmplex.write("prova.vtu", v);

        stats.stop_and_collect("vtu-write");

        int cell_num = 0;
        in.get("cell_num", cell_num);

        // Elem e;
        // space.elem(cell_num, e);
        // std::cout << e.measure() << std::endl;

        // V p0, p1, p2;
        // e.node(0, p0);
        // e.node(1, p1);
        // e.node(2, p2);

        // disp("-------------");
        // disp(p0);
        // disp("-------------");
        // disp(p1);
        // disp("-------------");
        // disp(p2);
        // disp("-------------");

        Quadrature q;
        auto grad = space.shape_grad(q);
        auto fun = space.shape(q);
        auto dx = space.differential(q);

        // DirichletBoundaryCondition<FunctionSpace> bc(space);

        PetscMatrix H;
        space.create_matrix(H);

        // disp(H);

        // Laplacian laplacian;

        {
            auto space_view = space.view_device();
            auto H_view = space.assembly_view_device(H);
            auto grad_view = grad.view_device();
            auto fun_view = fun.view_device();
            auto dx_view = dx.view_device();

            Device::parallel_for(space.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                Elem e;
                space_view.elem(i, e);

                auto g = grad_view.make(e);
                auto dx = dx_view.make(e);

                ElementMatrix el_mat;
                el_mat.set(0.0);

                for (PetscInt qp = 0; qp < Quadrature::NPoints; ++qp) {
                    for (PetscInt i = 0; i < Elem::NFunctions; ++i) {
                        for (PetscInt j = 0; j < Elem::NFunctions; ++j) {
                            el_mat(i, j) += inner(g(i, qp), g(j, qp)) * dx(qp);
                        }
                    }
                }

                // disp(g(0, 0));
                // disp("-------------");
                // disp(el_mat);
                space_view.add_matrix(e, el_mat, H_view);
            });
        }

        // disp(H);

        stats.stop_and_collect("assembly");
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(dmplex_test);
}  // namespace utopia

#endif  // UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)
