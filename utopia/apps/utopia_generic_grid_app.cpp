#include "utopia_StructuredGrid.hpp"
#include "utopia_ui.hpp"
#include "utopia_Views.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_Algorithms.hpp"

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
        constexpr double meas = g.measure();

        //COMPILE-TIME END
        assert(device::approxeq(meas, 1.0, 1e-8));

    }

    UTOPIA_REGISTER_APP(generic_grid_test);
}
