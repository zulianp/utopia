#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_Views.hpp"

#include <utility>

namespace utopia {

    class ViewTest {
    public:
        using Scalar = double;

        ViewTest()
        {}

        void run()
        {
            UTOPIA_RUN_TEST(array_view_test);
            UTOPIA_RUN_TEST(static_array_view_test);
            UTOPIA_RUN_TEST(vector_view_test);
            UTOPIA_RUN_TEST(matrix_view_test);
            UTOPIA_RUN_TEST(mv_view_test);
            UTOPIA_RUN_TEST(view_assign_test);
            UTOPIA_RUN_TEST(view_binary_test);
            UTOPIA_RUN_TEST(view_unary_test);
            UTOPIA_RUN_TEST(view_composite_test);
        }

        void array_view_test()
        {
            std::vector<Scalar> a_v(2, 1.0);
            ArrayView<Scalar> a(&a_v[0], a_v.size());
            device::axpy(4.0, a, a);
        }

        void static_array_view_test()
        {
            ArrayView<Scalar, 2> a;
            ArrayView<Scalar, 2, 2> b;

            device::fill(2.0, a);
            device::fill(3.0, b);

            device::axpy(4.0, a, a);

            device::axpy(4.0, b, b);
        }

        void vector_view_test()
        {
            Vector2<Scalar> a;
            a.set(1.0);

            Scalar dot_a = dot(a, a);
            utopia_test_assert(approxeq(2.0, dot_a));
        }

        void matrix_view_test()
        {
            StaticMatrix<Scalar, 2, 2> a;
            a.set(1.0);

            Scalar dot_a = dot(a, a);
            utopia_test_assert(approxeq(4.0, dot_a));
        }

        void mv_view_test()
        {
            StaticMatrix<Scalar, 3, 2> A;
            Vector2<Scalar> x;
            Vector3<Scalar> y;
            Vector3<Scalar> expected;
           
            A.set(1.0);
            x.set(2.0);

            expected.set(4.0);

            y = A * x;

            utopia_test_assert(approxeq(expected, y, 1e-10));
        }

        void view_assign_test()
        {
            using V = utopia::Vector3<Scalar>;
            V x;
            V y;
            x.set(0.0);
            y.set(1.0);

            DeviceAssign<V, V> va(x, y);
            va.apply();

            utopia_test_assert(approxeq(x, y));
        }

        void view_binary_test()
        {
            using V = utopia::Vector3<Scalar>;
            V x;
            V y;
            x.set(0.0);
            y.set(1.0);

            DeviceAssign<V, DeviceBinary<V, V, Plus>> va(x, DeviceBinary<V, V, Plus>(y, y));
            va.apply();

            y.scale(2.0);
            utopia_test_assert(approxeq(x, y));
        }

        void view_unary_test()
        {
            using V = utopia::Vector3<Scalar>;
            V x, y;
            x.set(2.0);
            y.set(4.0);

            x = pow2(x);
            x = abs(x);

            utopia_test_assert(approxeq(x, y));
        }

        void view_composite_test()
        {
            using V = utopia::Vector3<Scalar>;
            V x, y;

            x.set(2.0);
            y.set(6.0);

            x = pow2(x) + abs(x);

            utopia_test_assert(approxeq(x, y));
        }
    };

    void view()
    {
        ViewTest().run();
    }

    UTOPIA_REGISTER_TEST_FUNCTION(view);
}
