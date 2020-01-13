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
            UTOPIA_RUN_TEST(view_inv_test);
            UTOPIA_RUN_TEST(view_transpose_test);
            UTOPIA_RUN_TEST(view_trace_test);
            UTOPIA_RUN_TEST(view_norm_test);
            UTOPIA_RUN_TEST(view_eig_test);
            UTOPIA_RUN_TEST(view_diag_test);
            UTOPIA_RUN_TEST(choose_type);
            UTOPIA_RUN_TEST(size_test);
            UTOPIA_RUN_TEST(strain_test);
            UTOPIA_RUN_TEST(inner_test);
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
            StaticVector2<Scalar> a;
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
            StaticVector2<Scalar> x;
            StaticVector3<Scalar> y;
            StaticVector3<Scalar> expected;

            A.set(1.0);
            x.set(2.0);

            expected.set(4.0);

            y = A * x;

            utopia_test_assert(approxeq(expected, y, 1e-10));
        }

        void view_assign_test()
        {
            using V = utopia::StaticVector3<Scalar>;
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
            using V = utopia::StaticVector3<Scalar>;
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
            using V = utopia::StaticVector3<Scalar>;
            V x, y;
            x.set(2.0);
            y.set(4.0);

            x = pow2(x);
            x = abs(x);

            utopia_test_assert(approxeq(x, y));
        }

        void view_composite_test()
        {
            using V3 = utopia::StaticVector3<Scalar>;
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x3 = utopia::StaticMatrix<Scalar, 2, 3>;
            using Mat3x2 = utopia::StaticMatrix<Scalar, 3, 2>;
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;

            V3 x, y;

            x.set(2.0);
            y.set(6.0);

            x = pow2(x) + abs(x);

            utopia_test_assert(approxeq(x, y));

            Mat2x3 A;
            A.set(2.0);
            x.set(1.0);

            V2 w; w.set(1.0);
            V2 z = A * x + w;
            w.set(7.0);

            utopia_test_assert(approxeq(z, w));

            Mat3x2 A_t = transpose(A);

            Mat3x2 expected; expected.set(2.0);

            utopia_test_assert(approxeq(expected, A_t));

            A(0, 0) = 1.0;
            A(0, 1) = 2.0;
            A(0, 2) = 3.0;

            A(1, 0) = 3.0;
            A(1, 1) = 2.0;
            A(1, 2) = 1.0;

            Mat2x2 AAt = A * transpose(A);

            double det_A = det(AAt);
            utopia_test_assert(approxeq(96.0, det_A, 1e-8));
        }

        void view_transpose_test()
        {
            using Mat2x3 = utopia::StaticMatrix<Scalar, 2, 3>;
            using Mat3x2 = utopia::StaticMatrix<Scalar, 3, 2>;

            Mat2x3 A;

            A(0, 0) = 1.0;
            A(0, 1) = 2.0;
            A(0, 2) = 3.0;

            A(1, 0) = 3.0;
            A(1, 1) = 2.0;
            A(1, 2) = 1.0;


            Mat3x2 A_t = transpose(A);
            Mat3x2 A_t_expected;

            A_t_expected(0, 0) = 1.0;
            A_t_expected(1, 0) = 2.0;
            A_t_expected(2, 0) = 3.0;

            A_t_expected(0, 1) = 3.0;
            A_t_expected(1, 1) = 2.0;
            A_t_expected(2, 1) = 1.0;

            utopia_test_assert(approxeq(A_t, A_t_expected));
            utopia_test_assert(transpose(transpose(A)).is_alias(A));
        }

        void view_inv_test()
        {
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;

            Mat2x2 A;
            A(0, 0) = 1.0;
            A(0, 1) = 2.0;

            A(1, 0) = 3.0;
            A(1, 1) = 2.0;

            Mat2x2 A_inv = inv(A);
            Mat2x2 Id = A * A_inv;
            Mat2x2 expected_Id; expected_Id.set(0.0);

            expected_Id(0, 0) = 1.0;
            expected_Id(1, 1) = 1.0;

            utopia_test_assert(approxeq(Id, expected_Id, 1e-8));
            utopia_test_assert(inv(inv(A)).is_alias(A));
        }

        void view_trace_test()
        {
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;
            Mat2x2 A; A.set(2.0);
            utopia_test_assert(approxeq(trace(A), 4.0, 1e-8));
        }

        void view_norm_test()
        {
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x3 = utopia::StaticMatrix<Scalar, 2, 3>;

            V2 v; v.set(2.0);

            double n1 = norm1(v);
            double n2 = norm2(v);
            double n_infty = norm_infty(v);

            utopia_test_assert(approxeq(4.0, n1));
            utopia_test_assert(approxeq(std::sqrt(8.0), n2, 1e-8));
            utopia_test_assert(approxeq(2.0, n_infty));

            Mat2x3 m; m.set(2.0);

            n1 = norm1(m);
            n2 = norm2(m);
            n_infty = norm_infty(m);

            utopia_test_assert(approxeq(2*3*2.0, n1));
            utopia_test_assert(approxeq(std::sqrt(2*3*4.0), n2, 1e-8));
            utopia_test_assert(approxeq(2.0, n_infty));
        }

        void view_eig_test()
        {
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;
            Mat2x2 A; A.set(2.0);

            V2 e;

            eig(A, e);
            // disp(e);
        }

        void view_diag_test()
        {
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;
            Mat2x2 A; A.set(2.0);

            V2 d = diag(A);

            V2 expected; expected.set(2.0);

            utopia_test_assert(approxeq(d, expected));
        }

        void choose_type()
        {
            StaticMatrix<Scalar, 2, 2> A;

            auto a_t = transpose(A);
            using AT = decltype(a_t);


            static_assert(Traits<AT>::Order == 2, "must be 2nd order tensor");

            auto expr = transpose(A) + A;
            using E = decltype(expr);

            DeviceNumber<Scalar> num;
            using T = ChooseType<DeviceNumber<Scalar>, E, E>::Type;

            static_assert(Traits<T>::Order == 2, "must be 2nd order tensor");


            auto axA = 0.5 * A;
            using AXA = decltype(axA);
            static_assert(Traits<AXA>::Order == 2, "must be 2nd order tensor");
        }

        void size_test()
        {
            StaticMatrix<Scalar, 2, 2> A;
            rows( transpose(A) );
            utopia_test_assert((rows( transpose(A) + A )) == 2);
        }

        void strain_test()
        {
            StaticMatrix<Scalar, 2, 2> A, E;
            auto expr = 0.5 * (transpose(A) + A);
            E = expr;

            static_assert(Traits<decltype(expr)>::Order == 2, "must be 2nd order tensor");
        }

        void inner_test()
        {
            StaticMatrix<Scalar, 2, 2> A, E;
            auto expr = 0.5 * (transpose(A) + A);
            utopia_test_assert( rows(expr) == 2 );
        }
    };

    void view()
    {
        ViewTest().run();
    }

    UTOPIA_REGISTER_TEST_FUNCTION(view);
}
