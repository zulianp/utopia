
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include "utopia_Testing.hpp"

//include edsl components
#include "utopia_Core.hpp"

#include "utopia_kokkos_VectorView.hpp"
#include "utopia_kokkos_Traits.hpp"

#include <cmath>

namespace utopia {

    static void kokkos_vector_view()
    {
        using ViewType = Kokkos::View<double *>;

        ViewType kokkos_x("x", 10);
        VectorView<ViewType> x(kokkos_x);
        x.set(1.0);
        x.axpy(2., x);

        x += 0.5 * x;
        x *= 2.0;

        const double x_dot_x = dot(x, x);

        disp(x_dot_x);
    }

    static void kokkos_view()
    {
        UTOPIA_RUN_TEST(kokkos_vector_view);
    }

    UTOPIA_REGISTER_TEST_FUNCTION(kokkos_view);
}

#endif //WITH_TRILINOS
