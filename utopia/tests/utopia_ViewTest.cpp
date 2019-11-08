#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_ViewTraits.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_MatrixView.hpp"
#include "utopia_Algorithms.hpp"

#include <utility>

namespace utopia {

    class ViewTest {
    public:

        ViewTest()
        {}

        void run()
        {
            UTOPIA_RUN_TEST(array_view_test);
            UTOPIA_RUN_TEST(static_array_view_test);
            UTOPIA_RUN_TEST(vector_view_test);
            UTOPIA_RUN_TEST(matrix_view_test);
            UTOPIA_RUN_TEST(mv_view_test);
        }

        void array_view_test()
        {
            std::vector<double> a_v(2, 1.0);
            ArrayView<double> a(&a_v[0], a_v.size());
            device::axpy(4.0, a, a);
        }

        void static_array_view_test()
        {
            ArrayView<double, 2> a;
            ArrayView<double, 2, 2> b;

            device::fill(2.0, a);
            device::fill(3.0, b);

            device::axpy(4.0, a, a);

            device::axpy(4.0, b, b);
        }

        void vector_view_test()
        {
            VectorView< ArrayView<double, 2> > a;
            a.set(1.0);

            double dot_a = dot(a, a);
            utopia_test_assert(approxeq(2.0, dot_a));
        }

        void matrix_view_test()
        {
            MatrixView< ArrayView<double, 2, 2> > a;
            a.set(1.0);

            double dot_a = dot(a, a);
            utopia_test_assert(approxeq(4.0, dot_a));
        }

        void mv_view_test()
        {
            StaticMatrix<double, 3, 2> A;
            Vector2<double> x;
            Vector3<double> y;
            Vector3<double> expected;
           
            A.set(1.0);
            x.set(2.0);

            expected.set(4.0);

            y = A * x;

            utopia_test_assert(approxeq(expected, y, 1e-10));
        }
    };

    void view()
    {
        ViewTest().run();
    }

    UTOPIA_REGISTER_TEST_FUNCTION(view);
}
