#include "utopia_StructuredGrid.hpp"
#include "utopia_ui.hpp"
#include "utopia_Views.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_petsc_DMDA.hpp"

namespace utopia {
    using V = utopia::StaticVector<double, 2>;
    using I = utopia::ArrayView<int, 2>;
    using VA = utopia::ArrayView<double, 2>;
    template class StructuredGrid<V, I>;

    void generic_grid_test(Input &in)
    {
        //COMPILE-TIME BEGIN (everything is evaualted at compile time)
        constexpr StructuredGrid<V, I> g(
            {10, 10},
            {0, 0}, {10, 10},
            {0, 0}, {10, 10},
            VA({0.0, 0.0}),
            VA({1.0, 1.0})
        );

        static_assert(g.dim() == 2, "dim must be constexpr");
        static_assert(g.n_nodes() == 100, "n_nodes has to be computable at compile time");
        static_assert(g.n_elements() == 81, "wrong number of elements");
        constexpr double meas = g.measure();
        constexpr bool is_b = g.is_node_on_boundary(0);
        static_assert(is_b, "node 0 must be on boundary");
        // constexpr bool is_b = g.is_node_on_boundary_local_no_ghost(0, SideSet::left());

        //COMPILE-TIME END
        assert(device::approxeq(meas, 1.0, 1e-8));

    }

    UTOPIA_REGISTER_APP(generic_grid_test);

    void vec_grid_test(Input &in)
    {
        //run-time sizes real data
        std::vector<int>   zeros_v  = {0, 0};
        std::vector<int>   dims_v   = {10, 10};
        std::vector<float> box_min_v = {0.0f, 0.0f};
        std::vector<float> box_max_v = {1.0f, 1.0f};


        //views on data
        ArrayView<int>   zeros(&zeros_v[0], zeros_v.size());
        ArrayView<int>   dims(&dims_v[0], dims_v.size());
        ArrayView<float> box_min(&box_min_v[0], box_min_v.size());
        ArrayView<float> box_max(&box_max_v[0], box_max_v.size());

        StructuredGrid<VectorView<ArrayView<float>>, ArrayView<int>> g(
            dims,
            zeros, dims,
            zeros, dims,
            box_min,
            box_max
        );

        assert(g.dim() == 2);
        assert(g.n_nodes() == 100);
        assert(g.n_elements() == 81);
        double meas = g.measure();
        bool is_b = g.is_node_on_boundary(0);
        assert(is_b);

        assert(device::approxeq(meas, 1.0, 1e-8));

    }

    UTOPIA_REGISTER_APP(vec_grid_test);


    void dmda_test(Input &in)
    {
        PetscCommunicator comm;
        PetscDMDA<V, I> dmda(comm);

        dmda.read(in);

        PetscVector v;
        dmda.create_vector(v);

        dmda.write("prova.vtr", v);
    }

    UTOPIA_REGISTER_APP(dmda_test);
}
